# QC of the data
# 1. overall quality of the measurement
# 2. effects of spin before measurement
# Author: tao.xu
###############################################################################

setwd("D:/Users/tao.xu/Dropbox/Nurse project/")

data = read.csv("2013-04-10_Conc_Rui_Rabstein_Urin.csv")
data= data[-c(1:),]

## overall CV and within plate CV
index.ref = sapply(data$Sample.Identification, function(x) grep("Ref",x,fixed=T) )
index.ref = sapply(index.ref, function(x) length(x)!=0)

rst=NULL
for (i in measures){
	tmp = tapply(data[index.ref,i], INDEX = 		data$Plate.Bar.Code[index.ref], function(y) 		sd(y,na.rm=T)/mean(y,na.rm=T)
	)
	overall = sd(data[index.ref,i],na.rm=T)/mean(data[index.ref,i],na.rm=T)
	tmp = c(tmp,overall)
	rst = cbind(rst,tmp)
}
colnames(rst) = measures
#rst = rst[-1,2*(1:162)-1]
write.csv(t(rst), file = "CV of each metabolites across all paltes.csv")

## between plate differences of refsamples 
pdf("reference samples in different plates.pdf", width=15, height=15)
par(mfrow = c(5,5))
for(i in measures){
	data$m = data[,i]
	#plotmeans(m ~ Plate.Bar.Code, data, subset=index.ref)
	boxplot(m ~ Plate.Bar.Code,data,subset=index.ref, main = i)
}
dev.off()

## LOD
aboveLOD = function(data){
	index.zero = (data$Sample.Type=="Zero Sample")
	index.nurse = sapply(data$Sample.Identification, function(x) grep("SW",x,fixed=T) )
	index.nurse = sapply(index.nurse, function(x) length(x)!=0)
	
	rst.overLOD=apply(data[,measures],2, function(x) 	x>3*median(x[index.zero])
	)
	return(rst.overLOD)
}

rst=NULL; matrixLOD = data[,measures]
for(i in names(table(data$Plate.Bar.Code))[2:7]){
			subset = which(data$Plate.Bar.Code == i)
			matrixLOD[subset,]=aboveLOD(data[subset,])
}

		
aboveLOD = function(data){
	index.zero = (data$Sample.Type=="Zero Sample")
	index.nurse = sapply(data$Sample.Identification, function(x) grep("SW",x,fixed=T) )
	index.nurse = sapply(index.nurse, function(x) length(x)!=0)
	
	rst.overLOD=apply(data[,measures],2, function(x) 	sum(x[index.nurse]>3*median(x[index.zero]), na.rm = T)
	)
	return(rst.overLOD)
}

rst=NULL
for(i in names(table(data$Plate.Bar.Code))[2:7]){
	subset = which(data$Plate.Bar.Code == i)
	rst = cbind(rst,aboveLOD(data[subset,]))
}
colnames(rst) = names(table(data$Plate.Bar.Code))[2:7]
rst = data.frame(rst, overall = apply(rst,1,sum))
#for(i in 1:6){
#	rst[,i] = rst[,i]/table(data$Plate.Bar.Code[index.nurse])[i+1]
#}
rst[,7]=rst[,7]/429
#rst.overLOD = rst.overLOD[1:162)-1]
write.csv(rst, file = "overLOD_nurse.csv")


## data Normalization
tmp = apply(data[,measures], 2, function(x) 1000*x/data$Creatinine)
data[,measures] = tmp

##Spined and unspined samples
index.spin = sapply(data$Sample.Identification, function(x) grep("PU",x,fixed=T) )
index.spin = sapply(index.spin, function(x) length(x)!=0)

index.unspin = sapply(data$Sample.Identification, function(x) grep("PX",x,fixed=T) )
index.unspin = sapply(index.unspin, function(x) length(x)!=0)

measures = colnames(rst)
rst.difference = apply(data[,measures], 2, function(x) t.test(x[index.spin],x[index.unspin],paired=T))

rst = NULL
for(i in 1:length(rst.difference)){
	tmp = c(rst.difference[[i]]$estimate,
	rst.difference[[i]]$conf.int,
	rst.difference[[i]]$p.value)
	rst = rbind(rst,tmp)
}
rownames(rst) = measures
colnames(rst) = c("mean of difference","2.5%","97.5%","p-value")
write.csv(rst, "differences between spin and unspin.csv")

pdf("spin and unspined2.pdf", width = 15, height=15)
par(mfrow=c(5,5))
for(i in measures){
	plot(data[index.spin,i], type ="b", ylim = range(data[c(index.spin|index.unspin),i]), main = i)
	points(data[index.unspin,i], type ="b", col = "red")
}
dev.off()
#[2:8][c(-1,-9)]