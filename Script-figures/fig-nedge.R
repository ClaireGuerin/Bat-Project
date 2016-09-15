###############################################
###                         				###
### 		BAT IBM SENSORY JAMMING 		###
### 		EXPLORATORY STATISTICS			### 
###                         				###
###############################################

rm(list = ls())

library(car)
library(extrafont)
loadfonts()

###----------USER INPUT REQUIRED----------###
### Remember to create Scatterplots and Boxplots folders
resDir = "F:/Bats2016/Nedge/Nedge_variation_good/" # directory where the results are stored

varParam = "nedge"
dutyCycle = ""
# dutyCycle = "Unconstrained duty cycle" # uncomment accordingly
# dutyCycle = "5.4% duty cycle" # uncomment accordingly
# dutyCycle = " 15 duty cycle%" # uncomment accordingly

###----------###################----------###

startTime = Sys.time()

pcomb = c("nedge","tres","sim_durn","cornerpos","iid_on_axis","mov_angle","flight_speed",
	"call_durn","ipi","hear_tresh","source_level","alpha","speedsound")
pNames = c("Group size","Time Resolution (s)","Simulation Duration","Corner position",
	"IID on axis (m)","Movement angle","Flight speed","Call duration (s)",
	"Inter-pulse interval (s)","Hearing threshold (dB SPL)","Source level (dB SPL at 10 cm)",
	"Atmospheric absorption alpha","Speed of sound (m/s)")
varParIndex = which(pcomb == varParam)
Nindex = which(pcomb == "nedge")
varName = pNames[varParIndex]
indexNames = c("Echo-sound overlaps \n(mean number per echo per individual)","Echo-sounds total \noverlap fraction (mean %)",
	"% of echoes with overlaps", "% time of IPI free of sounds")

setwd(resDir)
repFolders = dir(resDir, pattern = "Res")
nRep = length(repFolders)
nSim = length(dir(repFolders[1]))

eelist = list()
eclist = list()
paramlist = list()
Nlist = list()

baseline = read.csv("baseline.csv", header=T, row.name=1, dec=".", sep=",")
varBase = which(colnames(baseline) == varParam)

#baseline = matrix(NA, nrow=902, ncol=10)
#colnames(baseline) = c(paste(c("meanNumOvPerEcho","meanTimeOvPerEch","X.echoOverlapped"),"C", sep=""),
#	paste(c("meanNumOvPerEcho","meanTimeOvPerEch","X.echoOverlapped"),"E", sep=""),pcomb[c(1,5,9,11)])

#baseline[,7] = rep(3,902)
#baseline[,8] = rep(1,902)
#baseline[,9] = rep(0.07,902)
#baseline[,10] = rep(120,902)
#summary(baseline)
#write.csv(baseline, file="baseline.csv")

color <- c("firebrick1","deepskyblue4")
color_transparent30 <- adjustcolor(color, alpha.f = 0.3) 
color_transparent20 <- adjustcolor(color, alpha.f = 0.2)
color_baseline <- c("darkred","cyan")

subplots = c("A","B","C")

for (i in 1:nSim){

###----------IMPORT INDICES----------###

	eelist[[i]] = matrix(NA, ncol = 7, nrow = 2)
	eclist[[i]] = matrix(NA, ncol = 7, nrow = 2)
	
	for (folder in repFolders){
		sim = dir(folder)[i]
		ee_file = paste(resDir,folder,"/",sim,"/","echo_echo_index.csv", sep="")
		ec_file = paste(resDir,folder,"/",sim,"/","echo_call_index.csv", sep="")
		param_file = paste(resDir,folder,"/",sim,"/","parameter_set.csv", sep="") 

		ee_ind = read.table(ee_file, header=T, sep=",", dec=".", row.names = 1)
		ec_ind = read.table(ec_file, header=T, sep=",", dec=".", row.names = 1)
		parameter = as.numeric(as.character(read.table(param_file, header=F, sep="\n", dec=".")[varParIndex,]))
		N = as.numeric(as.character(read.table(param_file, header=F, sep="\n", dec=".")[Nindex,])) ^ 2

		colnames(eelist[[i]]) = colnames(ee_ind)
		colnames(eclist[[i]]) = colnames(ec_ind)
		
		eelist[[i]] = rbind(eelist[[i]], ee_ind)
		eclist[[i]] = rbind(eclist[[i]], ec_ind)
		paramlist[[i]] = parameter
		Nlist[[i]] = N
	}

###----------End of INDICES IMPORT----------###	

###----------BOXPLOTS----------###

	ec_ids = as.factor(as.vector(eclist[[i]][1])[,1])
	ee_ids = as.factor(as.vector(eelist[[i]][1])[,1])

	imgPath1 = paste("Boxplots/",sim,".pdf", sep="")
	pdf(imgPath1, family="CM Roman", width=7, height=7)
	layout(matrix(c(1,1,2,3,4,5), ncol=2, byrow=TRUE), widths=c(1,1,1), heights=c(1,4,4))

	par(mai=c(0,0,0,0))
	plot(c(0,2), c(0,2), pch=20, col="white", axes=F, ylab="", xlab="")
	text(1,1,labels=paste(varName,"=",as.character(paramlist[[i]])), cex=3)

	par(mai=c(0.6,0.8,0.42,0.42))
	for (j in 1:3){
		
		indexCol = c(2,4,6)[j]
		index = indexNames[j]
		ec_index = as.vector(eclist[[i]][indexCol])[,1]
		ee_index = as.vector(eelist[[i]][indexCol])[,1]
		n = Nlist[[i]] - 1
		
		if (j==1){
			indlim = n 
		}else{
			indlim = max(100, ec_index, ee_index, na.rm=T)
		}

		Boxplot(ec_index~ec_ids, main=subplots[j], col.main="gray30", col=color_transparent30[1], border=color[1], cex.axis=1.5, ylab = index, xlab = "Individual ID", ylim=c(0,indlim), cex.lab=1.5, col.lab="gray30", id.method="none") 
		Boxplot(ee_index~ee_ids, col=color_transparent30[2], border=color[2], add=T, ylab = "", xlab = "", axes=F, id.method="none")
	}
	par(mai=c(0,0,0,0))
	plot.new()
	legend(x="center", ncol=1, legend=c("Echo-call overlaps","Echo-echo overlaps"), text.col = "gray30",
      	fill=color_transparent30, border=color, title=expression(italic("Type of overlapping sound")), cex=2, bty="n")
	
	dev.off()

###----------End of BOXPLOTS----------###

}

###----------SCATTERPLOTS----------###

imgPath2 = paste("Scatterplots/",varParam,substr(dutyCycle,1,3),".pdf", sep="")
pdf(imgPath2, family="CM Roman", width=7, height=7)
layout(matrix(c(1,1,2,3,4,5), ncol=2, byrow=TRUE), widths=c(1,1,1), heights=c(1,4,4))

par(mai=c(0,0,0,0))
plot(c(0,2), c(0,2), pch=20, col="white", axes=F, ylab="", xlab="")
text(1,1,labels=paste(varName,"\n",dutyCycle), cex=3)

par(mai=c(0.6,0.8,0.42,0.42))	
for (j in 1:3){

	indexCol = c(2,4,6)[j]
	index = indexNames[j]
	ec_index_all = c()
	ec_index_means = c()
	ee_index_means = c()
	ee_index_all = c()
	param_all = c()
	param_for_means = c()
	n_all = c()
	
	for (i in 1:nSim){

		ec_index = as.vector(eclist[[i]][indexCol])[,1]
		ee_index = as.vector(eelist[[i]][indexCol])[,1]
		len_ind = length(ec_index)
		param = rep(Nlist[[i]], len_ind)
		n = Nlist[[i]] - 1

		ec_index_all = append(ec_index_all, ec_index)
		ec_index_means = append(ec_index_means, mean(ec_index, na.rm=T))
		ee_index_all = append(ee_index_all, ee_index)
		ee_index_means = append(ee_index_means, mean(ee_index, na.rm=T))
		param_all = append(param_all, param)
		param_for_means = append(param_for_means,Nlist[[i]])
		n_all = append(n_all, n)
		ec_means = cbind(ec_index_means,param_for_means)
		ec_means = ec_means[order(param_for_means), ]	
		ee_means = cbind(ee_index_means,param_for_means)
		ee_means = ee_means[order(param_for_means), ]	
	}

	if (j==1){
		maxy = max(n_all) 
	}else{
		maxy = max(100, ec_index_all, ee_index_all, na.rm=T)
	}	
	jitterStrength = 1

	#baseline[,j] = ec_index_all
	#baseline[,j+3] = ee_index_all

	plot(ec_index_all~jitter(param_all, jitterStrength), main=subplots[j], col.main="gray30", cex.axis=1.5, col=color_transparent30[1], ylab = index, 
		xlab = varName, ylim=c(0,maxy), cex.lab=1.5, col.lab="gray30", xaxt="n", pch=21, bg=color_transparent20[1]) 
	axis(1, at = signif(param_for_means, digits=3), las=1)
	points(ee_index_all~jitter(param_all, jitterStrength), cex.axis=1.5, col=color_transparent30[2], pch=21, bg=color_transparent20[2])
	points(baseline[,j]~jitter(baseline[,varBase]^2,jitterStrength), cex.axis=1.5, col=color_baseline[1], pch=8)
	points(baseline[,j+3]~jitter(baseline[,varBase]^2,jitterStrength), cex.axis=1.5, col=color_baseline[2], pch=8)
	points(ec_index_means~param_for_means, data=ec_means, col=color[1], type="c", lwd=1.5)
	points(ee_index_means~param_for_means, data=ee_means, col=color[2], type="c", lwd=1.5)
}
par(mai=c(0,0,0,0))
plot.new()
legend(x="center", ncol=1, 
	legend=c("Echo-call overlaps","Echo-call overall mean","Echo-call for baseline parameters",
	"Echo-echo overlaps","Echo-echo overall mean","Echo-echo for baseline parameters"),
	col=c(color[1],color_baseline[1],color_baseline[1],color[2],color_baseline[2], color_baseline[2]), text.col = "gray30",
	title=expression(italic("Type of overlapping sound")), pch=c(21,NA,8,21,NA,8), lty=c(NA,1,NA,NA,1,NA), pt.cex=3, 
	cex=1.5, pt.bg=c(color_transparent30[1],NA,NA,color_transparent30[2],NA,NA), bty="n", lwd=c(NA,1.5,NA,NA,1.5,NA))		
dev.off()

###----------End of SCATTERPLOTS----------###

print(Sys.time() - startTime)