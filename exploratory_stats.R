###############################################
###                         				###
### 		BAT IBM SENSORY JAMMING 		###
### 		EXPLORATORY STATISTICS			### 
###                         				###
###############################################

rm(list = ls())

library(Rcmdr)

resDir = "D:/Bat_Project/Res/"
setwd(resDir)
repFolders = dir(resDir, pattern = "Res")
nRep = length(repFolders)
nSim = length(dir(repFolders[1]))

eelist = list()
eclist = list()

for (i in 1:nSim){

	eelist[[i]] = matrix(NA, ncol = 7, nrow = 2)
	colnames(eelist[[i]]) = c("ID","meanNumOvPerEcho","sdNumOv","meanTimeOvPerEch","sdTimeOv","%echoOverlapped","%timeFreeIPI")
	eclist[[i]] = matrix(NA, ncol = 7, nrow = 2)
	colnames(eclist[[i]]) = c("ID","meanNumOvPerEcho","sdNumOv","meanTimeOvPerEch","sdTimeOv","%echoOverlapped","%timeFreeIPI")
	
	for (folder in repFolders){
		sim = dir(folder)[i]
		ee_file = paste(resDir,folder,"/",sim,"/","echo_echo_index.csv", sep="")
		ec_file = paste(resDir,folder,"/",sim,"/","echo_call_index.csv", sep="")

		ee_ind = read.table(ee_file, header=T, sep=",", dec=".", row.names = 1)
		ec_ind = read.table(ec_file, header=T, sep=",", dec=".", row.names = 1)

		eelist[[i]] = rbind(eelist[[i]], ee_ind)
		eclist[[i]] = rbind(eclist[[i]], ec_ind)
	}
	
	ec_ids = as.factor(as.vector(eclist[[i]][1])[,1])
	ee_ids = as.factor(as.vector(eelist[[i]][1])[,1])

	imgPath = paste("Exploratory/",sim,".pdf", sep="")
	pdf(imgPath, width=7, height=7)
	par(mfrow = c(2,2))

	indexNames = c("Mean number of sound overlaps per echo","Mean time of sound overlap per echo", "Number of overlapped echoes (%)", "Time of IPI free of sounds (%)")
	
	for (j in c(2,4,6,7)){
		
		index = indexNames[j]
		ec_index = as.vector(eclist[[i]][j])[,1]
		ee_index = as.vector(eelist[[i]][j])[,1]
		indlim = max(ec_index, ee_index, na.rm=T)

		Boxplot(ec_index~ec_ids, medcol="coral", cex=3, ylab = index, xlab = "ID", ylim=c(0,indlim)) 
		Boxplot(ee_index~ee_ids, medcol="mediumseagreen", cex=3, ylab = index, xlab = "ID", add=T)
	}
	
	dev.off()	
}