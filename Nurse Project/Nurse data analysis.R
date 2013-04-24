# 1. Data visualization
# 2. 
# Author: tao.xu
###############################################################################

samples = read.csv("2013-01-31 Helmholtz.csv")
samples = samples[1:429,]

index.person = match(samples$Proben_ID[which(samples$SW_Nr == "SW1030")], data$Sample.Identification)

pdf("personal change over time SW1032_unnormalized.pdf", width =20, height=15)
par(mfrow = c(5,5))
for(i in measures){
	plot(data[index.person, i], type = "b", main = i, pch = c(21,19)[matrixLOD[index.person,i]+1])
}
dev.off()

plotNurse = function(data, samples, id){
	index.person = match(samples$Proben_ID[which(samples$SW_Nr == id)], data$Sample.Identification)
	pdf(paste("personal change over time ",id,".pdf", sep = "",collapse=""), width =20, height=15)
	par(mfrow = c(5,5))
	for(i in measures){
		subset = as.logical(matrixLOD[index.person,i])
		if(sum(subset)!=0)	{
			plot(data[index.person[subset], i], type = "b", main = i, col = c("red", "blue")[samples$Schichtdienst[which(samples$SW_Nr == id)][subset]])
			points(data[index.person[subset], i], pch = 22, col = c("black")[samples$Morgenurin[which(samples$SW_Nr == id)][subset]])
		}
		else plot(0, main = i)
	}
	dev.off()
}

plotNurse = function(data, samples, id){
	index.person = match(samples$Proben_ID[which(samples$SW_Nr == id)], data$Sample.Identification)
	pdf(paste("personal change over time ",id,".pdf", sep = "",collapse=""), width =20, height=15)
	par(mfrow = c(5,5))
	for(i in measures){
		subset = as.logical(matrixLOD[index.person,i])
		if(sum(subset)!=0)	{
			plot(data[index.person[subset], i]~data[index.person[subset], "Probennahme_Uhr"], 
					main = i, 
					pch = c(19, 21)[samples$Schichtdienst[which(samples$SW_Nr == id)][subset]])
			#points(data[index.person[subset], i], pch = 22, col = c("black")[samples$Morgenurin[which(samples$SW_Nr == id)][subset]])
		}
		else plot(0, main = i)
	}
	dev.off()
}

plotNurse(data.merged,samples,"SW1037")

for(i in names(table(samples$SW_Nr))){
	plotNurse(data.merged,samples,i)
}


pdf(paste("Change over time in all samples",id,".pdf", sep = "",collapse=""), width =20, height=15)
	par(mfrow = c(5,5))
	for(i in measures){
		subset = as.logical(matrixLOD[,i])
		if(sum(subset)!=0)	{
			plot(data.merged[subset, i]~data.merged[subset, "Probennahme_Uhr"], #log = 'y',
					main = i, 
					pch = c(19, 21)[data.merged$Schichtdienst[subset]])
			#points(data[index.person[subset], i], pch = 22, col = c("black")[samples$Morgenurin[which(samples$SW_Nr == id)][subset]])
		}
		else plot(0, main = i)
	}
	dev.off()


## differences between night shift and day shift work
require(nlme)

tmp = merge(data, samples, by.x = "Sample.Identification", by.y = "Proben_ID")
data.merged = tmp
matrixLOD = matrixLOD[-304,]

rst=NULL
for(i in measures){
	subset = as.logical(matrixLOD[83:510,i])
	if(sum(subset)>3&
			table(data.merged$Schichtdienst[subset])[1]!=0&
			table(data.merged$Schichtdienst[subset])[2]!=0
	){
		data.merged$m = data.merged[,i]
		model = lme(m ~ Schichtdienst, data.merged[subset,], random = ~1|SW_Nr)
		rst = rbind(rst, summary(model)$tTable[2,])
	}
	else rst = rbind(rst,rep(NA,5))
}
rownames(rst) = measures
rst = data.frame(rst, fdr = p.adjust(rst$p.value, method = "BH"), bonf = p.adjust(rst$p.value, method = "bonf"))

