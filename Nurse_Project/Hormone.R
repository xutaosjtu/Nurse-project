hormone = read.csv("data/Hormone data_sul.csv" )
colnames(hormone)[10:12] = c("Melatonin", "Cortisol", "Estradiol")

hormone = hormone[order(hormone$SW_Nr, hormone$Probennahme_Dat, hormone$Proben_Nr),]
#hormone$Probennahme_Uhr = hormone$Probennahme_Uhr*24
hormone$Probennahme_Dat=as.character(hormone$Probennahme_Dat)

hormone$Melatonin = scale(log(hormone$Melatonin))
hormone$Cortisol = scale(log(hormone$Cortisol))
hormone$Estradiol = scale(log(hormone$Estradiol))

hormone$time = paste(hormone$Probennahme_Dat, hormone$Probennahme_Uhr)
hormone$time = strptime(hormone$time, "%Y.%m.%d %H:%M")

metabo = data.merged
metabo[, valid_measures] = scale(log(metabo[, valid_measures]))

for(p in levels(hormone$SW_Nr)){
  
  metabo_person = subset(metabo, SW_Nr==p&Schichtdienst=='Tagschicht')
  hormone_person = subset(hormone, SW_Nr==p&Schichtdienst=='Tagschicht')
    
  yrange = c(-3,3)
  
  plot(metabo_person$time, metabo_person$C0, type = 'b', col="red", ylim = yrange, pch = 1)
  points(hormone_person$time, hormone_person$Cortisol, type = 'b', col = "red", pch = 2)
  points(hormone_person$time, hormone_person$Melatonin, type = 'b', col = "red", pch = 3)
  points(hormone_person$time, hormone_person$Estradiol, type = 'b', col = "red", pch = 4)
  
  metabo_person = subset(data.merged, SW_Nr==p&Schichtdienst=='Nachtschicht')
  hormone_person = subset(hormone, SW_Nr==p&Schichtdienst=='Nachtschicht')
  
  
  yrange = c(0, max(max(metabo_person$C0,  na.rm=T), max(hormone_person$Melatonin, na.rm=T)))
  
  plot(metabo_person$time, metabo_person$C0, type = 'b', col="blue", ylim = yrange, pch = 1)
  points(hormone_person$time, hormone_person$Cortisol, type = 'b', col = "blue", pch = 2, lty=2)
  points(hormone_person$time, hormone_person$Melatonin, type = 'b', col = "blue", pch = 3, lty=3)
  points(hormone_person$time, hormone_person$Estradiol, type = 'b', col = "blue", pch = 4, lty=4)
}



##ToDo: plot function for the correlation between metabolite and hormone levels
for(p in unique(data.merged$SW_Nr)){  
  for(shift in c('Tagschicht', 'Nachtschicht')){
    metabo_person = subset(metabo, SW_Nr==p&Schichtdienst==shift)
    hormone_person = subset(hormone, SW_Nr==p&Schichtdienst==shift)
    
    y1 = aggregate.2(metabo_person$C0, by=list(metabo_person$hour_cut,metabo_person$Probennahme_Dat))
    y2 = aggregate.2(hormone_person$Cortisol, by=list(hormone_person$hour_cut, hormone_person$Probennahme_Dat))
    
    if(shift == "Tagschicht"){
      col = "red"
    }
    else {col="blue"}
    plot(y1$time, y1$x, type="b", ylim = range(y1$x, y2$V1), col=col, pch = 2)
    points(y2$time, y2$V1, type="b", col=col, pch = 3)
  }
  
}


aggregate.2 = function(data, by){
  y = aggregate(data, by, mean, na.rm=T)
  y = as.data.frame(y)
  y[,1] = as.character(y[,1])
  y[,1] = sapply(y[,1], function(x) return(unlist(strsplit(x, split=" ", fixed=T))[2]))
  y$time = paste(y$Group.2, y$Group.1)
  y$time = strptime(y$time, "%Y.%m.%d %H:%M")
  return(y)
}



#############################################################################
###   Analysis:
###   1. Categorize the sampling time, calculate the aggregated hormone levels and metabolite concentrations at different time intervals
###   2. Correlation analysis between metabolite concentrations and hormone levels.
###   3. Analyze the metabolite secretion rate in categorized time intervals.
###   
#############################################################################

### categorized metabolite concentration and hormone levels
### time intervals were set at 0:00:00~8:00:00, 8:00:00~16:00:00, 16:00:00~ 24:00:00
hormone$hour = strptime(hormone$Probennahme_Uhr, "%H:%M")
hormone$hour_cut = cut(hormone$hour, breaks = strptime(c("0:00:00","8:00:00", "16:00:00", "23:59:59"), "%H:%M:%S"))

data.merged$hour = strptime(data.merged$Probennahme_Uhr, "%H:%M")
data.merged$hour_cut = cut(data.merged$hour, breaks = strptime(c("0:00:00","8:00:00", "16:00:00", "23:59:59"), "%H:%M:%S"))
metabo = data.merged
metabo[, valid_measures] = scale(log(metabo[, valid_measures]))


## aggregate the metabolite concentration according to the time intervales
metabo_aggre = aggregate(metabo[,valid_measures], by=list(metabo$SW_Nr, metabo$hour_cut,metabo$Probennahme_Dat, metabo$Schichtdienst), mean, na.rm=T)
colnames(metabo_aggre)[1:4]=c("SW_Nr","hour", "date", "shift")

hormone_aggre = aggregate(hormone[, c("Melatonin","Cortisol","Estradiol")], by=list(hormone$SW_Nr, hormone$hour_cut, hormone$Probennahme_Dat, hormone$Schichtdienst), mean, na.rm=T)
colnames(hormone_aggre)[1:4]=c("SW_Nr","hour", "date", "shift")

## matched metabolite measurements and hormone levels
match_metabo_hormone = merge(metabo_aggre, hormone_aggre)

match_metabo_hormone = match_metabo_hormone[order(match_metabo_hormone$SW_Nr, match_metabo_hormone$date, match_metabo_hormone$hour),]

## heatmap of correlation between metabolites and hormorns
require(corrplot)
pdf("correlation between metabolite with hormones.pdf", width = 20, height=20)
col1 = colorRampPalette(c("#053061","#2166AC","#4393C3","#92C5DE", "#D1E5F0", "#FFFFFF", "#FDDBC7","#F4A582","#D6604D", "#B2182B",  "#67001F"))  
corrplot(cor(match_metabo_hormone[,-c(1:4)], use="pair"), tl.col = "black")
dev.off()

## using linear mixed effect model to estimate the association between hormone and metabolites
require(nlme)
rst = NULL
for(i in valid_measures){
  match_metabo_hormone$m = match_metabo_hormone[,i]
  model = lme(m~ Melatonin + Cortisol + Estradiol, data = match_metabo_hormone, random = ~1|SW_Nr, na.action = na.omit)
  tmp = summary(model)$tTable
  rst = rbind(rst, c(tmp[1,], tmp[2, ], tmp[3,]))
}
rownames(rst) = valid_measures
write.csv(rst, "association between metabolites and hormones_lme.csv")

## 3. Analyze the metabolite secretion rate in categorized time intervals.
metabo_aggre = merge(metabo_aggre, samples.addition, by.x = "SW_Nr", by.y = "P_ID")
metabo_aggre = metabo_aggre[order(metabo_aggre$SW_Nr, metabo_aggre$date, metabo_aggre$hour),]
levels(metabo_aggre$Kontrolle)=c("case", "control")

##GEE
require(gee)
rst = NULL
for(i in valid_measures){
  metabo_aggre$m = metabo_aggre[,i]
    model = gee(m ~ as.factor(Kontrolle)+ Alter + BMI + as.factor(AR_Rauch_zurzt)+as.factor(SD), #  
                id = SW_Nr, 
                data = metabo_aggre, 
                subset = shift=="Tagschicht"& SW_Nr!="SW1041", #&data.merged$Alter>=45
                na.action=na.omit, 
                corstr = "exchangeable"
    )
    rst = rbind(rst, summary(model)$coef[2,])
}
rownames(rst) = valid_measures
rst = data.frame(rst, p.value = 2*pnorm(-abs(rst[,5])))
rst = data.frame(rst, fdr = p.adjust(rst$p.value, method = "BH"), bonf = p.adjust(rst$p.value, method = "bonf"))
rst$Estimate = -rst$Estimate
write.csv(rst, "categorized time_Chronic effect of night shift work_GEE_daywork_age_BMI_smoking_disease_exclude diab.csv")

require(nlme)
rst=NULL
for(i in valid_measures){
  metabo_aggre$m = metabo_aggre[,i]
  model = lme(m ~ shift + Alter + BMI + as.factor(AR_Rauch_zurzt) + as.factor(SD),
              metabo_aggre,
                subset = Kontrolle=="case" & SW_Nr!="SW1041",
                random = ~ 1|SW_Nr,
                na.action=na.omit
  )
  rst = rbind(rst, summary(model)$tTable[2,])
}
rownames(rst) = valid_measures
rst = data.frame(rst)
rst$Value = -rst$Value
rst = data.frame(rst, fdr = p.adjust(rst$p.value, method = "BH"), bonf = p.adjust(rst$p.value, method = "bonf"))
write.csv(rst, file = "categorized time_Short term effect of night shift_mixed model_age_BMI_smoking_disease.csv")


plot(metabo_aggre$C0[which(metabo_aggre$SW_Nr=="SW1030")], type = "b", col=c("blue", "red")[metabo_aggre$shift[which(metabo_aggre$SW_Nr=="SW1030")]])

plot(match_metabo_hormone$C0[which(match_metabo_hormone$SW_Nr=="SW1030")], type = "b", col=c("blue", "red")[match_metabo_hormone$shift[which(match_metabo_hormone$SW_Nr=="SW1030")]])

plot(C0~hour, match_metabo_hormone, subset=(shift=="Tagschicht"&Kontrolle=="case"))
plot(C0~hour, match_metabo_hormone, subset=shift=="Nachtschicht"&Kontrolle=="case", add = T, col = "grey")

plot(C0~hour, match_metabo_hormone, subset=(shift=="Tagschicht"&Kontrolle=="case"))
plot(C0~hour, match_metabo_hormone, subset=shift=="Tagschicht"&Kontrolle=="control", add = T, col = "grey")

require(ggplot2)
require(grid)

pdf("metabolite concentration at different time period of the day.pdf", width = 30, height = 20)
k=1;nrow=1;ncol=1;newp=F
pushViewport(viewport(layout = grid.layout(3, 3)))
for(i in valid_measures){
  match_metabo_hormone$m = match_metabo_hormone[,i]
  p = ggplot(match_metabo_hormone, aes(hour, m))
  p = p + geom_boxplot(aes(fill=interaction(factor(Kontrolle), factor(shift))))+ggtitle(i)
  
  print(p, vp = viewport(layout.pos.row = nrow, layout.pos.col = ncol))
  if(ncol==3 & nrow==3) {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(3, 3)))
  }
  if(k<9) {
    if(k%%3==0) nrow=nrow+1
    k = k+1
    ncol = k-3*(nrow-1)
  }
  else {
    k=1
    nrow = 1
    ncol = 1
  }
    
}
dev.off()

p = ggplot(match_metabo_hormone, aes(x = Cortisol, y = C0, group = shift, col =factor(shift)))
p + geom_point() + geom_smooth()# + geom_line(aes(y = Cortisol), col="red")


pdf("metabolite correlation with Cortisol.pdf", width = 20, height = 20)
k=1;nrow=1;ncol=1;newp=F
pushViewport(viewport(layout = grid.layout(3, 3)))
for(i in valid_measures){
  match_metabo_hormone$m = match_metabo_hormone[,i]
  
  ##change the code here if plot some thing else
  p = ggplot(match_metabo_hormone, aes(x = Cortisol, y = m))
  p = p + geom_point() + geom_smooth(method = "lm")+ggtitle(i)
  ##
  
  print(p, vp = viewport(layout.pos.row = nrow, layout.pos.col = ncol))
  if(ncol==3 & nrow==3) {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(3, 3)))
  }
  if(k<9) {
    if(k%%3==0) nrow=nrow+1
    k = k+1
    ncol = k-3*(nrow-1)
  }
  else {
    k=1
    nrow = 1
    ncol = 1
  }
  
}
dev.off()
