# 1. Data visualization: Show metabolite changes overtime
# 2. Find differences between day and night shift work
# Author: tao.xu
###############################################################################

index.person = match(samples$Proben_ID[which(samples$SW_Nr == "SW1030")], data$Sample.Identification)

## Data visualization
pdf("personal change over time SW1032_unnormalized.pdf", width =20, height=15)#Change over time of unnormalized metabolite concentrations in one sample
par(mfrow = c(5,5))
for(i in measures){
	plot(data[index.person, i], type = "b", main = i, pch = c(21,19)[matrixLOD[index.person,i]+1])
}
dev.off()

plotNurse = function(data, matrixLOD, id){
	## plot metabolite concentrations overtime for each nurse  
	## input: data: metabolite measuresments of M meatbolites in N sample (N*M+x), with x other variables
	##		  matrixLOD: matrix indicating the measurements above LOD (N*M+y), with y other variables
	## 		  id: ID number of the nurse
	index.person = which(matrixLOD$SW_Nr==id)
	pdf(paste("personal change over time ",id,".pdf", sep = "",collapse=""), width =20, height=15)
	par(mfrow = c(5,5))
	for(i in measures){
		subset = as.logical(matrixLOD[index.person,i])
		if(sum(subset)!=0 & sum(!is.na(data[index.person[subset],i]))!=0){
			tmp = data[index.person[subset],]
			tmp = tmp[order(tmp$Proben_Nr),]
			tmp$m = tmp[,i]
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
				#points(data[index.person[subset], i], pch = 22, col = c("black")[samples$Morgenurin[which(samples$SW_Nr == id)][subset]])
				nl=nl+1
			}	
		}
		else plot(0, main = i)
	}
	dev.off()
}

plotNurse(data.merged, matrixLOD,"SW1030")

for(i in names(table(samples$SW_Nr))){
	plotNurse(data.merged,matrixLOD,i)
}


pdf("Change over time in all samples.pdf", width =20, height=15)# Change over time in all samples
	par(mfrow = c(5,5))
	for(i in measures){
		subset = as.logical(matrixLOD[,i])
		if(sum(subset)!=0)	{
			plot(data.merged[subset, i]~data.merged[subset, "Probennahme_Uhr"], 
					#log = 'y',
					main = i, 
					pch = c(19, 21)[data.merged$Schichtdienst[subset]])
			#points(data[index.person[subset], i], pch = 22, col = c("black")[samples$Morgenurin[which(samples$SW_Nr == id)][subset]])
		}
		else plot(0, main = i)
	}
dev.off()


## Find differences between night shift and day shift work
# Linear mixed effect model
require(nlme)
rst=NULL
for(i in measures){
	subset = as.logical(matrixLOD[,i])
	if(sum(subset)>3& i!="Creatinine"&
			table(data.merged$Schichtdienst[subset])[1]!=0&
			table(data.merged$Schichtdienst[subset])[2]!=0
	){
		data.merged$m = as.numeric(as.character(data.merged[,i]))
		model = lme(m ~ Schichtdienst, data.merged[which(subset),], random = ~1|SW_Nr, na.action=na.exclude)
		rst = rbind(rst, summary(model)$tTable[2,])
	}
	else rst = rbind(rst,rep(NA,5))
}
rownames(rst) = measures
rst = data.frame(rst)
rst = data.frame(rst, fdr = p.adjust(rst$p.value, method = "BH"), bonf = p.adjust(rst$p.value, method = "bonf"))


## Find chronic effects of night shift and day shift work
# Linear mixed effect model
data.merged$group = rep(1,dim(data.merged)[1])
data.merged$group[which(data.merged$SW_Nr=="SW1039"|data.merged$SW_Nr=="SW1041")]=0
data.merged$group=as.factor(data.merged$group)
require(nlme)
rst=NULL
for(i in measures){
  subset = as.logical(matrixLOD[,i])
  data.merged$m = data.merged[,i]
  if(sum(subset)>3& i!="Creatinine"&table(data.merged$group[subset])[1]!=0&
       table(data.merged$group[subset])[2]!=0){
    #model = lme(m ~ group, data.merged, random = ~1|SW_Nr,na.action=na.exclude)
    #rst = rbind(rst, summary(model)$tTable[2,])
    model = lm(m ~ group, data.merged, na.action=na.exclude)
    rst = rbind(rst, summary(model)$coef[2,])
  }
  else rst = rbind(rst,rep(NA,5))
}
rownames(rst) = measures
rst = data.frame(rst)
rst = data.frame(rst, fdr = p.adjust(rst$p.value, method = "BH"), bonf = p.adjust(rst$p.value, method = "bonf"))

##GEE
rst = NULL
for(i in measures[-c(16,39, 81:140)]){
  subset = as.logical(matrixLOD[,i])
  data.merged$m = scale(log(data.merged[,i]))
  if(sum(subset)>100& i!="Creatinine"&table(data.merged$group[subset])[1]!=0&
       table(data.merged$group[subset])[2]!=0){
    model = gee(group ~ m, id = SW_Nr, data = data.merged, subset = which(data.merged$Schichtdienst=="Tagschicht"), na.action=na.omit, corstr = "exchangeable", family = binomial)
    rst = rbind(rst, summary(model)$coef[2,])
  }
  else rst = rbind(rst,rep(NA,5))
}
rownames(rst) = measures[-c(16,39, 81:140)]
rst = data.frame(rst, pvalue = 2*pnorm(-abs(rst[,5])))
write.csv(rst, "Chronic effect of night shift work_GEE_daywork.csv")