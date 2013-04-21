# QC of the data
# 1. overall quality of the measurement
# 2. effects of spin before measurement
# Author: tao.xu
###############################################################################

setwd("D:/Users/tao.xu/Dropbox/Nurse project/")

data = read.csv("2013-04-10_Conc_Rui_Rabstein_Urin.csv")
data= data[-c(1:),]

index.ref = sapply(data$Sample.Identification, function(x) grep("Ref",x,fixed=T) )
index.ref = sapply(index.ref, function(x) length(x)!=0)

rst=NULL
for (i in 19:342){
	tmp = tapply(data[index.ref,i], INDEX = 		data$Plate.Bar.Code[index.ref], function(y) 		sd(y,na.rm=T)/mean(y,na.rm=T)
	)
	overall = sd(data[index.ref,i],na.rm=T)/mean(data[index.ref,i],na.rm=T)
	tmp = c(tmp,overall)
	rst = cbind(rst,tmp)
}
colnames(rst) = colnames(data)[19:342]
rst = rst[-1,2*(1:162)-1]
write.csv(t(rst), file = "CV of each metabolites across all paltes.csv")


index.zero = (data$Sample.Type=="Zero Sample")
##calculation of LOD
#LOD=...
#

index.nurse = sapply(data$Sample.Identification, function(x) grep("SW",x,fixed=T) )
index.nurse = sapply(index.nurse, function(x) length(x)!=0)

rst.overLOD=apply(data[,19:342],2, function(x) sum(x[index.nurse]>x[1]))
rst.overLOD=rst.overLOD/429
rst.overLOD = rst.overLOD[2*(1:162)-1]
write.csv(rst.overLOD, file = "overLOD_nurse.csv")


##
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


