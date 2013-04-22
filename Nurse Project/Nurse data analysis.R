# 1. Data visualization
# 2. 
# Author: tao.xu
###############################################################################

samples = read.csv("2013-01-31 Helmholtz.csv")
samples = samples[1:429,]

index.person = which(data$Sample.Identification %in% samples$Proben_ID[which(samples$SW_Nr == "SW1032")])

pdf("personal change over time SW1032_unnormalized.pdf", width =20, height=15)
par(mfrow = c(5,5))
for(i in measures){
	plot(data[index.person, i], type = "b", main = i, pch = c(21,19)[matrixLOD[index.person,i]+1])
}
dev.off()
