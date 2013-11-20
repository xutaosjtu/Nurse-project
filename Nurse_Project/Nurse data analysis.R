# 1. Data visualization: Show metabolite changes overtime
# 2. Find differences between day and night shift work
# Author: tao.xu
###############################################################################

require(nlme)
require(gee)

pdf("correlation between metabolites.pdf", width = 10, height = 10)
col3 <- colorRampPalette(c("red", "white", "blue"))  
corrplot(cor(data.merged[,valid_measures], use = "pair"),col=col3(100), tl.col="black")
dev.off()

data.merged$Schichtdienst =factor(data.merged$Schichtdienst, levels = c("Tagschicht", "Nachtschicht"))
data.merged = data.merged[order(data.merged$SW_Nr, data.merged$Probennahme_Dat, data.merged$Proben_Nr),]

#data.merged$Probennahme_Uhr = sapply(data.merged$Probennahme_Uhr, function(x) return(x*24))
data.merged$Probennahme_Dat = as.character(data.merged$Probennahme_Dat)
data.merged$time = paste(data.merged$Probennahme_Dat, data.merged$Probennahme_Uhr)
data.merged$time = strptime(data.merged$time, "%Y.%m.%d %H:%M")

### calculate the time interval between two samples
data.merged$time_interval = 0
for(i in 1:(nrow(data.merged)-1)){
  if(data.merged$Probennahme_Uhr[i]<data.merged$Probennahme_Uhr[i+1]){
    data.merged$time_interval[i+1] = data.merged$Probennahme_Uhr[i+1]-data.merged$Probennahme_Uhr[i]
  }
  else data.merged$time_interval[i+1] = data.merged$Probennahme_Uhr[i+1]+24-data.merged$Probennahme_Uhr[i]
}


pdf(paste("Nurse metabolite concentration cumulative all.pdf"), width = 60, height=160)
par(mfcol = c(58,25))
for(p in levels(data.merged$SW_Nr)){
  sample_person = which(data.merged$SW_Nr==p)
  tmp_person = data.merged[sample_person,]
  for(m in valid_measures){
    nl = 0
    sums=tapply(tmp_person[,m]*tmp_person$time_interval, INDEX = tmp_person$Probennahme_Dat, sum, na.rm=T)
    for(day in unique(tmp_person$Probennahme_Dat)){
      sample_day = which(data.merged$Probennahme_Dat[sample_person]==day)
      tmp = tmp_person[sample_day,]
      
      tmp[which(is.na(tmp[,m])),m] = mean(tmp[,m], na.rm = T)
      tmp[,m] = tmp[,m]*tmp$time_interval
      
      tmp.rst = cumsum(tmp[,m])
      if(nl==0){
        plot(tmp.rst~ tmp$Probennahme_Uhr[which(tmp$Probennahme_Dat==day)], type = "b", xlim = c(0,24), ylim = range(0,max(sums)*1.5), main = m, col = tmp$Schichtdienst[1])
      }
      else points(tmp.rst~ tmp$Probennahme_Uhr[which(tmp$Probennahme_Dat==day)], type = "b", xlim = c(0,24), col = tmp$Schichtdienst[1])
      nl = nl+1
    }
  }
}
dev.off()


# index.person = match(samples$Proben_ID[which(samples$SW_Nr == "SW1030")], data$Sample.Identification)
# 
# ## Data visualization
# pdf("personal change over time SW1032_unnormalized.pdf", width =20, height=15)#Change over time of unnormalized metabolite concentrations in one sample
# par(mfrow = c(5,5))
# for(i in measures){
# 	plot(data[index.person, i], type = "b", main = i, pch = c(21,19)[matrixLOD[index.person,i]+1])
# }
# dev.off()

plotNurse = function(data, id){
	## plot metabolite concentrations overtime for each nurse  
	## input: data: metabolite measuresments of M meatbolites in N sample (N*M+x), with x other variables
	##		  matrixLOD: matrix indicating the measurements above LOD (N*M+y), with y other variables
	## 		  id: ID number of the nurse
	tmp = subset(data, SW_Nr==id)
	pdf(paste("personal change over time ",id,".pdf", sep = "",collapse=""), width =20, height=15)
	par(mfrow = c(5,5))
	for(i in valid_measures){
	  tmp$m = tmp[,i]
		#subset = as.logical(matrixLOD[index.person,i])
		if(sum(!is.na(tmp$m))!=0){
			#tmp = data[index.person[subset],]
			#tmp = tmp[order(tmp$Proben_Nr),]
			Levels = unique(tmp[,"Probennahme_Dat"])
			nl=0;
			for(l in Levels){
				indexL=which(tmp$Probennahme_Dat==l)
				if(nl == 0){
					plot(m~Probennahme_Uhr,
							tmp[indexL,],
							#subset = ,
							main = i, type = "b", 
							pch = c(19,21)[tmp$Schichtdienst[indexL]],
							col = c("red","blue")[tmp$Schichtdienst[indexL]],
							xlim = c(0,23), ylim = range(tmp$m, na.rm = T)
					)
				}
				else{
					points(m~Probennahme_Uhr,
							tmp,
							subset = which(tmp$Probennahme_Dat==l),
							main = i, type = "b", 
							pch = c(19,21)[tmp$Schichtdienst[indexL]],
							col = c("red","blue")[tmp$Schichtdienst[indexL]]
					)
				}
				nl=nl+1
			}	
		}
		else plot(0, main = i)
	}
	dev.off()
}

plotNurse(data.merged, "SW1036")

for(i in names(table(samples$SW_Nr))){
	plotNurse(data.merged,i)
}

subset = data.merged$SW_Nr=="SW1040"
plot(log(data.merged[subset, i]), pch = c(19, 21)[data.merged$Schichtdienst[subset]],
     col = c("red", "blue")[data.merged$Schichtdienst[subset]], type = 'b')


pdf("Change over time in all samples.pdf", width =20, height=15)# Change over time in all samples
	par(mfrow = c(5,5))
	for(i in measures){
		subset = !is.na(data.merged[,i])
		if(sum(subset)!=0)	{
			plot(log(data.merged[subset, i])~data.merged[subset, "Probennahme_Uhr"], 
					#log = 'y',
					main = i, 
					pch = c(19, 21)[data.merged$Schichtdienst[subset]],
           col = c("red", "blue")[data.merged$Schichtdienst[subset]]
           )
			#points(data[index.person[subset], i], pch = 22, col = c("black")[samples$Morgenurin[which(samples$SW_Nr == id)][subset]])
		}
		else plot(0, main = i)
	}
dev.off()


## Compare the cumulative metabolite secretion in Urin
tmp = data.merged

data_columative = NULL
for(i in valid_measures){
  tmp[,i] = unlist(tapply(tmp[,i], INDEX=tmp$Probennahme_Dat, function(x) {x[which(is.na(x))]=mean(x, na.rm=T);return(x)}))
  if(is.null(data_columative)) data_columative = aggregate(tmp[,i]*tmp$time_interval, by=list(tmp$SW_Nr, tmp$Probennahme_Dat), sum, na.rm=T)
  else  data_columative = cbind(data_columative, aggregate(tmp[,i]*tmp$time_interval, by=list(tmp$SW_Nr, tmp$Probennahme_Dat), sum, na.rm=T)[,3])
}
colnames(data_columative) = c("SW_Nr", "Probennahme_Dat", valid_measures)

data_columative = merge(data_columative, 
          samples[!duplicated(samples[, c("SW_Nr", "Probennahme_Dat")]),]
          )
data_columative = merge(data_columative, samples.addition, by.x = "SW_Nr", by.y = "P_ID")
data_columative = data_columative[order(data_columative$SW_Nr, data_columative$Probennahme_Dat),]

rst = NULL
for(i in valid_measures){
  data_columative$m = data_columative[,i]
  model = gee(m ~ as.factor(Kontrolle)+ Alter + BMI + as.factor(AR_Rauch_zurzt), # + as.factor(batch) +as.factor(SD)
              id = SW_Nr, 
              data = data_columative, 
              subset = Schichtdienst=="Tagschicht"& SW_Nr!="SW1041", #&Alter>=45
              na.action=na.omit, 
              corstr = "exchangeable"
  )
  rst = rbind(rst, summary(model)$coef[2,])
}
rownames(rst) = valid_measures
rst = data.frame(rst, p.value = 2*pnorm(-abs(rst[,5])))
rst = data.frame(rst, fdr = p.adjust(rst$p.value, method = "BH"), bonf = p.adjust(rst$p.value, method = "bonf"))


#   for(p in levels(tmp$SW_Nr)){
#     sample_person = which(data.merged$SW_Nr==p)
#     tmp_person = data.merged[sample_person,]
#     for(m in valid_measures){
#       nl = 0
#       sums=tapply(tmp_person[,m]*tmp_person$time_interval, INDEX = tmp_person$Probennahme_Dat, sum, na.rm=T)
#     }
#   }
}

samples[!duplicated(samples[,c("Schichtdienst","SW_Nr")]), ]

## Find chronic effects of night shift and day shift work
data.merged$group = rep(1,nrow(data.merged))
data.merged$group[which(data.merged$Kontrolle=="K")]=0
data.merged$group=as.factor(data.merged$group)

batch2 = names(table(data$Plate.Bar.Code))[7:12]
data.merged$batch = rep(1, nrow(data.merged))
data.merged$batch[which(data.merged$Plate.Bar.Code %in% batch2)] = 2
data.merged$batch = as.factor(data.merged$batch)
rm(batch2)

# Linear mixed effect model
# rst=NULL
# for(i in valid_measures){
#   subset = as.logical(matrixLOD[,i])
#   data.merged$m = data.merged[,i]
#   if(sum(subset)>100& i!="Creatinine"&table(data.merged$group[subset])[1]!=0&
#        table(data.merged$group[subset])[2]!=0){
#     model = lme(m ~ group, data.merged, random = ~1|SW_Nr,na.action=na.exclude)
#     rst = rbind(rst, summary(model)$tTable[2,])
#     #model = lm(m ~ group, data.merged, na.action=na.exclude)
#     #rst = rbind(rst, summary(model)$coef[2,])
#   }
#   else rst = rbind(rst,rep(NA,5))
# }
# rownames(rst) = valid_measures
# rst = data.frame(rst)
# rst = data.frame(rst, fdr = p.adjust(rst$p.value, method = "BH"), bonf = p.adjust(rst$p.value, method = "bonf"))
# write.csv("")

##GEE
rst = NULL
for(i in valid_measures){
  subset = !is.na(data.merged[,i])
  data.merged$m = log(data.merged[,i])
  if(sum(subset)>100& i!="Creatinine" & table(data.merged$group[subset])[1]!=0&
       table(data.merged$group[subset])[2]!=0){
    model = gee(m ~ as.factor(group)+ Alter + BMI + as.factor(AR_Rauch_zurzt) + as.factor(batch) +as.factor(SD), #  
                id = SW_Nr, 
                data = data.merged, 
                subset = Schichtdienst=="Tagschicht"& SW_Nr!="SW1041", #&data.merged$Alter>=45
                na.action=na.omit, 
               # family = binomial,
                corstr = "exchangeable"
                )
    rst = rbind(rst, summary(model)$coef[2,])
  }
  else rst = rbind(rst,rep(NA,5))
}
rownames(rst) = valid_measures
rst = data.frame(rst, p.value = 2*pnorm(-abs(rst[,5])))
rst = data.frame(rst, fdr = p.adjust(rst$p.value, method = "BH"), bonf = p.adjust(rst$p.value, method = "bonf"))
write.csv(rst, "Chronic effect of night shift work_GEE_daywork_age_BMI_smoking_disease_batch_exclude diab.csv")

pdf("metabolite concentration between cases and controls.pdf", width =20, height=15)
par(mfrow = c(5,10))
for(i in measures){
  if(sum(!is.na(data.merged[,i])!=0)){
    plot(data.merged$group,data.merged[,i],log="y", main = paste(i, "case vs control"))
    plot(data.merged$Schichtdienst,data.merged[,i],log="y", main = paste(i, "day vs night"), col = "grey")
  }
  else {
    plot(0, main =paste(i, "case vs control"))
    plot(0, main =paste(i, "day vs night"))
  }
}
dev.off()

## Find differences between night shift and day shift work
# Linear mixed effect model
rst=NULL
for(i in valid_measures){
  subset = !is.na(data.merged[,i])
  if(sum(subset)>10& i!="Creatinine"&
       table(data.merged$Schichtdienst[subset])[1]!=0&
       table(data.merged$Schichtdienst[subset])[2]!=0
  ){
    data.merged$m = scale(log(data.merged[,i]))
    model = lme(m ~ Schichtdienst + Alter + BMI + AR_Rauch_zurzt+ as.factor(SD)+as.factor(batch),
                data.merged,
                subset = which(data.merged$group == 1),
                random = ~ 1|SW_Nr,
                na.action=na.exclude
    )
    rst = rbind(rst, summary(model)$tTable[2,])
  }
  else rst = rbind(rst,rep(NA,5))
}
rownames(rst) = valid_measures
rst = data.frame(rst)
rst$Value = -rst$Value
rst = data.frame(rst, fdr = p.adjust(rst$p.value, method = "BH"), bonf = p.adjust(rst$p.value, method = "bonf"))
write.csv(rst, file = "Short term effect of night shift_mixed model_age_BMI_smoking_disease_batch.csv")

##gee
rst=NULL
for(i in valid_measures){
  subset = !is.na(data.merged[,i])
  if(sum(subset)>10& i!="Creatinine"&
       table(data.merged$Schichtdienst[subset])[1]!=0&
       table(data.merged$Schichtdienst[subset])[2]!=0
  ){
    data.merged$m = scale(log(data.merged[,i]))
    model = gee( m ~ Schichtdienst + Alter+ BMI + AR_Rauch_zurzt + as.factor(SD) + Plate.Bar.Code
                 , id = SW_Nr
                 , data = data.merged
                 , subset = which(data.merged$group == 1)
                 , na.action=na.omit 
                 , corstr = "exchangeable"
    )
    rst = rbind(rst, summary(model)$coef[2,])
  }
  else rst = rbind(rst,rep(NA,5))
}
rownames(rst) = valid_measures
rst = data.frame(rst, pvalue = 2*pnorm(-abs(rst[,5])))
write.csv(rst, file = "Short term effect of night shift_GEE_age.csv")
