###############################################
###                         				###
### 		BAT IBM SENSORY JAMMING 		###
### 		EXPLORATORY STATISTICS			### 
###                         				###
###############################################

rm(list = ls())

library(Rcmdr)

###----------USER INPUT REQUIRED----------###

resDir = "D:/Bat_Project/Res/" # directory where the results are stored

varParam = "ipi"
# dutycycle = "Unconstrained" # uncomment accordingly
# dutycycle = "5.4%" # uncomment accordingly
# dutycycle = "15%" # uncomment accordingly

###----------###################----------###

startTime = Sys.time()

pcomb = c("nedge","tres","sim_durn","cornerpos","iid_on_axis","mov_angle","flight_speed",
	"call_durn","ipi","hear_tresh","source_level","alpha","speedsound")
varParIndex = which(pcomb == varParam)

setwd(resDir)
repFolders = dir(resDir, pattern = "Res")
nRep = length(repFolders)
nSim = length(dir(repFolders[1]))

eelist = list()
eclist = list()
paramlist = list()

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
		
		colnames(eelist[[i]]) = colnames(ee_ind)
		colnames(eclist[[i]]) = colnames(ec_ind)
		
		eelist[[i]] = rbind(eelist[[i]], ee_ind)
		eclist[[i]] = rbind(eclist[[i]], ec_ind)
		paramlist[[i]] = parameter
	}

###----------End of INDICES IMPORT----------###	

###----------BOXPLOTS----------###

	ec_ids = as.factor(as.vector(eclist[[i]][1])[,1])
	ee_ids = as.factor(as.vector(eelist[[i]][1])[,1])

	imgPath1 = paste("Boxplots/",sim,".pdf", sep="")
	pdf(imgPath1, width=7, height=7)
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

###----------End of BOXPLOTS----------###

}

###----------SCATTERPLOTS----------###

imgPath2 = paste("Scatterplots/",sim,".pdf", sep="")
pdf(imgPath2, width=7, height=7)
par(mfrow = c(2,2))
	
for (j in 1:4){

	indexCol = c(2,4,6,7)[j]
	index = indexNames[j]
	ec_index_all = c()
	ee_index_all = c()
	param_all = c()
	
	for (i in 1:nSim){

		ec_index = as.vector(eclist[[i]][indexCol])[,1]
		ee_index = as.vector(eelist[[i]][indexCol])[,1]
		len_ind = length(ec_index)
		param = rep(paramlist[[i]], len_ind)

		ec_index_all = append(ec_index_all, ec_index)
		ee_index_all = append(ee_index_all, ee_index)
		param_all = append(param_all, param)		
	}

	maxy = max(ec_index_all, ee_index_all, na.rm=T)
	jitterStrength = 0.1*(max(param_all)-min(param_all))
	
	plot(ec_index_all~jitter(param_all, jitterStrength), col="coral", pch=20, ylab = index, xlab = varParam, ylim=c(0,maxy)) 
	points(ee_index_all~jitter(param_all, jitterStrength), col="mediumseagreen", pch=20)
}
		
dev.off()

###----------End of SCATTERPLOTS----------###

print(Sys.time() - startTime)

layout(matrix(c(1,2,3,3), ncol=2, byrow=TRUE), heights=c(4, 1))

par(mai=rep(0.5, 4))
plot(1:3,4:6,main="plot 1")
plot(1:3,4:6,main="plot 2")

par(mai=c(0,0,0,0))
plot.new()
legend(x="center", ncol=3,legend=c("0-1 km","1-5 km","outside barrier"),
       fill=c("green","orange","red"), title="Fetch")