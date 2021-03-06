# QC of the data
# 1. Overall quality of the measurement
# 2. Data normalization
# 3. Effects of sample spinning
# 4. Batch effects after normalization
# Author: tao.xu
###############################################################################

setwd("D:/Users/tao.xu/Dropbox/Nurse project/")

#data = read.csv("data/2013-04-10_Conc_Rui_Rabstein_Urin_normalisiert.csv")
#data = read.csv("data/original data/2013-04-10_Conc_Rui_Rabstein_Urin.csv")
data = read.csv("data/Urin data_all samples.csv")
#data= data[-c(1:),]
#samples = read.csv("data/2013-01-31 Helmholtz.csv")
#samples = samples[1:429,]
samples = read.csv("data/Sample map_all sampel.csv")
samples.addition = read.csv("data/sample_additional infor.csv")
samples = merge(samples, samples.addition, by.x = "SW_Nr", by.y = "P_ID")
rm(samples.addition)

measures = colnames(data)[(9:170)*2+1]

###################################################
# Overall quality of the measurement
# 1. Overall CV and within plate CV
# 2. Between plate differences of refsamples
# 3. Measurements above LOD
###################################################
index.ref = sapply(data$Sample.Identification, function(x) grep("Ref",x,fixed=T) )
index.ref = sapply(index.ref, function(x) length(x)!=0)


##over all CV and within plate CV 
rst=NULL
for (i in measures){
	within = tapply(data[index.ref,i], 
                  INDEX = data$Plate.Bar.Code[index.ref], 
                  function(y) sd(y,na.rm=T)/mean(y,na.rm=T)
	)
	overall = sd(data[index.ref,i],na.rm=T)/mean(data[index.ref,i],na.rm=T)
	rsti = c(within,overall)
	rst = cbind(rst,rsti)
}
colnames(rst) = measures
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

## Measurements above LOD
aboveLOD = function(data){
	# measurement above LOD in one plate
	# input: meatbolite measurements of samples and zero samples in one plate
	# output: Logical indicator of above LOD
	index.zero = (data$Sample.Type=="Zero Sample")
	index.nurse = sapply(data$Sample.Identification, function(x) grep("SW",x,fixed=T) )
	index.nurse = sapply(index.nurse, function(x) length(x)!=0)
	print(sum(index.nurse))
	#rst.overLOD=apply(data[,measures],2, function(x) {x=rep(3*median(x[index.zero]),length(x)); return(x)})
	rst.overLOD=apply(data[,measures],2, function(x) x>3*median(x[index.zero], na.rm = T))
	return(rst.overLOD[index.nurse, ])
}

matrixLOD = data[,measures] # indicator matrix of measurements above LOD
for(i in names(table(data$Plate.Bar.Code))){
	plate = which(data$Plate.Bar.Code == i)
	matrixLOD[plate,]=aboveLOD(data[plate,])
}
matrixLOD = data.frame(matrixLOD, 
		Sample.Identification=data$Sample.Identification)
matrixLOD = merge(matrixLOD,samples, by.x = "Sample.Identification", by.y = "Proben_ID")

rst=NULL # measurements above LOD in each plate and all plates
for(i in names(table(data$Plate.Bar.Code))){
	subset = which(data$Plate.Bar.Code == i)
	tmp = aboveLOD(data[subset,])
  print(nrow(tmp))
	rst = cbind(rst,apply(tmp , 2, sum, na.rm = T))
}
colnames(rst) = names(table(data$Plate.Bar.Code))
rst = data.frame(rst, overall = apply(rst,1,sum))
rst$overall=rst$overall/nrow(samples)
write.csv(rst, file = "overLOD_nurse.csv")

## data Normalization
data.merged = merge(data,samples,by.x="Sample.Identification", by.y="Proben_ID")

normalize<-function(data, measures){
	tmp = apply(data[,measures], 2, function(x) 1000*x/data$Creatinine)
	data[,measures] = tmp
	return(data)
}

data.merged = data.merged[which(matrixLOD$Creatinine==1), ]# exclude abnormal creatinine data
matrixLOD = matrixLOD[which(matrixLOD$Creatinine==1),]# exclude abnormal creatinine data
data.merged = normalize(data.merged,measures)

for(i in measures){
	data.merged[which(matrixLOD[,i]==0),i]=NA
}

## outliers
index=apply(data.merged[,measures],2, function(x) which(abs(x)>mean(x,na.rm=T)+5*sd(x,na.rm=T)|abs(x)<mean(x,na.rm=T)-5*sd(x,na.rm=T))) 
for(i in measures){
	data.merged[index[[i]],i]=NA
}



## Effects of sample spinning
index.spin = sapply(data$Sample.Identification, function(x) grep("PU",x,fixed=T) )
index.spin = sapply(index.spin, function(x) length(x)!=0)

index.unspin = sapply(data$Sample.Identification, function(x) grep("PX",x,fixed=T) )
index.unspin = sapply(index.unspin, function(x) length(x)!=0)

data = normalize(data,measures)

measures = colnames(rst)
rst.difference = apply(data[,measures], 2, function(x) t.test(x[index.spin],x[index.unspin],paired=T)) ## t-test
rst = NULL #differences between spin and unspined samples
for(i in 1:length(rst.difference)){
	tmp = c(rst.difference[[i]]$estimate,
	rst.difference[[i]]$conf.int,
	rst.difference[[i]]$p.value)
	rst = rbind(rst,tmp)
}
rownames(rst) = measures
colnames(rst) = c("mean of difference","2.5%","97.5%","p-value")
write.csv(rst, "differences between spin and unspin.csv")

pdf("spin and unspined.pdf", width = 15, height=15) # scatter plot to show trends of meatbolites 
par(mfrow=c(5,5))
for(i in measures){
	plot(data[index.spin,i][2:8], type ="b", ylim = range(data[c(index.spin|index.unspin),i][c(-1,-9)]), main = i, col = "red")
	points(data[index.unspin,i][2:8], type ="b", col = "blue")
}
dev.off()
#[2:8][c(-1,-9)]


## Batch effects after normalization
pdf("Metabolite levels in different plates_log.pdf", width = 25, height = 20)
par(mfrow=c(5,5))
for(i in measures){
	subset = as.logical(matrixLOD[83:510,i])
	if(sum(subset)!=0){
		plot(data.merged[subset,i]~data.merged$Plate.Bar.Code[subset], log="y", main=i)
	}
	else plot(0)
}
dev.off()

pdf("Metabolite levels of reference samples, internal standard and zero samples in different plates_log.pdf", width = 25, height = 20)
par(mfrow=c(5,5))
for(i in measures){
	#subset = as.logical(matrixLOD[83:510,i])
		plot(data[index.ref,i]~data$Plate.Bar.Code[index.ref], main=i, ylim = c(min(data[index.zero,i]),max(data[index.ref,i])))
		lines(data$Plate.Bar.Code[index.standard],data[index.standard,i],  main=i)
		boxplot(data[index.zero,i]~data$Plate.Bar.Code[index.zero],add=T, col = "blue")
}
dev.off()
