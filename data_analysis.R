###############################################
###                         				###
### 		BAT IBM SENSORY JAMMING 		###
### 		DESCRIPTIVE STATISTICS			### 
###                         				###
###############################################

### Created on 22/06/2016 at 15:15
### Author: Claire

rm(list = ls())
# clear current workspace

library(IRanges)

startTime = Sys.time()

resDir = "C:/Users/tbeleyur/Documents/Bat-Project/IPI_15DC/"
setwd(resDir)
resFiles = dir(resDir, pattern="Res")

for (resfile in resFiles){

resNum = resfile
combFolders = dir(resNum)

for (comb in combFolders){

newDir = paste(resDir,resNum,"/",comb, sep="")
# directory where the specific results that you want
# to extract are stored. Should be a folder containing 3 sub-folders:
# Calling, Moving, Hearing, with the corresponding data

#---------- PARAMETERS OF THE SIMULATION ----------#
# This section extracts from the parameter combination used in
# a the simulation, that is to be analysed.

param = read.csv(paste(newDir,"/","parameter_set.csv",sep=""), header = F)

SIMDURATION <- as.numeric(as.character(param[3,1])) # duration of the simulation iterations
TIMERESOLUTION <- as.numeric(as.character(param[2,1])) # time resolution, in seconds per iteration
FLIGHTSPEED <- as.numeric(as.character(param[7,1])) # speed of flight of the bats, in m/s
VSOUND <- as.numeric(as.character(param[13,1])) # m/s
TETA <- as.numeric(as.character(param[6,1])) # bat movement angle, in radians
SOURCELEVEL <- as.numeric(as.character(param[11,1])) # source level of the bat's call in dB SPL at 10 cm
HEARINGTHRESHOLD <- as.numeric(as.character(param[10,1])) # hearing threshold of the bat in dB SPL
ALPHA <- as.numeric(as.character(param[12,1])) # sound absorbtion rate at particular frequency in db/m  

#---------- End of SIMULATION PARAMETERS ----------#

#---------- DATA IMPORT & RESHAPE ----------#
# This section imports the data from the simulation 
# ran under python with the sensory_jamming.py program.
# It also reshapes the data into a R-friendly data format.

### CALLING ###

C.path = paste(newDir,"/","Calling/",sep="") # set path to the Calling folder.
C.fileNames = dir(C.path, pattern =".txt") 
# names of all files in the folder in alphabetical order, 
# i.e.: file names are ordered in IDs ascending order.

cal = as.data.frame(matrix(NA, ncol = length(C.fileNames), nrow = SIMDURATION))
# empty data frame of dimensions simulation duration * number of bats.

for (i in 1:length(C.fileNames)){ 
# for every file in the Calling folder
	cal[,i] = read.table(paste(C.path,C.fileNames[i], sep=""), header = F, sep = "\t", dec = ".")
	# import data into the data frame
	colnames(cal)[i] = paste("Bat",i-1, sep="")
	# calling information: 1 column per bat x, named Batx
}

### MOVING ###

M.path = paste(newDir,"/","Moving/",sep="") # set path to the Moving folder.
M.fileNames = dir(M.path, pattern =".txt") 
# names of all files in the folder in alphabetical order, 
# i.e.: file names are ordered in IDs ascending order, with
# x-coordinates data first and y-coordinates second.

mov = as.data.frame(matrix(NA, ncol = length(M.fileNames), nrow = SIMDURATION))
# empty data frame of dimensions simulation duration * (2*number of bats).
id_M = ceiling(1:length(M.fileNames)/2) # bat identification per column 
colnames(mov) = paste("Bat",id_M-1,c(".X",".Y"), sep = "")
# calling information: 2 column per bat x, named Batx.X, Batx.Y 

for (i in id_M){
# for every file in the Moving folder
	mov[,1+2*(i-1)] = read.table(paste(M.path,M.fileNames[1+2*(i-1)], sep=""), header = F, sep = "\t", dec = ".")
	# import x-coordinate data into the data frame
	mov[,2+2*(i-1)] = read.table(paste(M.path,M.fileNames[2+2*(i-1)], sep=""), header = F, sep = "\t", dec = ".")
	# import y-coordinate data into the data frame
}

### HEARING ###

H.path = paste(newDir,"/","Hearing/",sep="") # set path to the Hearing folder.
H.fileNames = dir(H.path, pattern =".txt")
# names of all files in the folder in alphabetical order, 
# i.e.: file names are ordered in IDs ascending order, with
# c-data first, i-data second and t-data last.
# See README.txt for information on c, i & t-data.

for (i in seq(1,length(H.fileNames),3)){
# for every bat in the simulation 
   	
	ownID = (ceiling(1:length(H.fileNames)/3)-1)[i] # ID of the bat
	tcall = read.table(paste(H.path,H.fileNames[i], sep=""), header = F, sep = "\t", dec = ".") # import c-data
	idbat = read.table(paste(H.path,H.fileNames[i+1], sep=""), header = F, sep = "\t", dec = ".") # import i-data
	tmstp = read.table(paste(H.path,H.fileNames[i+2], sep=""), header = F, sep = "\t", dec = ".") # import t-data

	temp = as.data.frame(matrix(NA, ncol = 4, nrow = dim(tcall)[1]))
	# temporary empty data frame to store data, of dimensions number of heared items * 3
	colnames(temp) = c("ID_R","Time_R","ID_E","Time_E")
	# 3 columns corresponding to the 3 data types c, i and t.
	temp$ID_R = ownID
	temp$Time_R = tmstp[,1] # transfer t-data into temp
	temp$ID_E = idbat[,1] # transfer i-data into temp
	temp$Time_E = tcall[,1] # transfer c-data into temp
	own_temp = temp[which(temp$ID_E == ownID),] 
	# extract sounds heard, which were originally emitted by the bat itself
	other_temp = temp[-which(temp$ID_E == ownID),]
	# extract sounds heard, which were originally emitted by other bats
	assign(paste("calls.self",ownID,sep=""),own_temp) 
	# name of object: calls.selfx for bat x
	assign(paste("calls.others",ownID,sep=""),other_temp) 
	# name of object: calls.othersx for bat x
}

#---------- End of DATA IMPORT & RESHAPE ----------#

#---------- ECHOES CALCULATION ----------#
# Calculate the time of arrival of primary echoes.
# Calculations of time of reception of the echo are only made 
# for echoes from calls that were produced by the bat itself.

### R.DIST FUN ###
# Function calculating the time at which the travelling sound,
# which bounces-off after reaching another individual than the 
# emittor (i.e. the bat who called in the first place), "collides" with
# each other of the moving bats.

R.dist = function(expl.t){ 
# R.dist depends on a unique variable: time.
# expl.t stands for exploratory time. expl.t is typically a vector of 
# time delays, for which possible bat-sound "collisions" are tested.
# expl.t is typically comprised between 0 and the time that was needed
# for the call to travel from its source to the bat on which
# it has bounced-off, thus producing an echo.

	new.x = init.x + FLIGHTSPEED * expl.t * TIMERESOLUTION * cos(TETA)
	# new position of the moving bat on the x-axis, after delta t = expl.t time
	new.y = init.y + FLIGHTSPEED * expl.t * TIMERESOLUTION * sin(TETA)
	# new position of the moving bat on the y-axis, after delta t = expl.t time
	R.t.sq = (new.x - x.source) ** 2 + (new.y - y.source) ** 2
	# R.t.sq stands for R(t)². Distance between the moving bat, which
	# originally emitted the sound, and the position where the echo 
	# bounced-off, otherwise called "source" as in source of the echo.
	# R(t)² slowly decreases, then increases with time.
	r.t.sq = (VSOUND * expl.t*TIMERESOLUTION) ** 2
	# r.t.sq stands for r(t)². Radius of the sound as it propagates,
	# i.e. distance between the source of the echo and its current
	# position in space. r(t)² increases with time.
	difrnce<-abs(R.t.sq - r.t.sq	)
	collision<- which(difrnce==min(difrnce))
	# a bat-sound collision occurs as soon as R(t)² <= r(t)²
	# Exact time of impact
	return(collision)
}

### HEAR.ECHO FUNCTION ###
# Function to calculate the max hearing distance for an echo 
# Inputs : 
# preCDist: distance sound travels from calling bat to target bat before impact/collision (metres)
# postCDist: distance sound travels from target bat - now as an echo - towards the bat which had called (metres)
# SL : source level of the call  (dB SPL at 10cm)
# hearTh: hearing threshold of a bat (dB SPL)

# output:
# Boolean T/F value on whether a bat can hear the echo or not - based on the given hearing threshold

Hear.echo = function(preCDist,postCDist, SL, hearTh){
	batRadius = 0.1 # bat simplified as a sphere of 10 cm radius
	batTS = 20*log10(batRadius/2) # target strength of the bat 
	SL_1m = SL - abs(20*log10(0.1/1)) # Source level calculated at 1 m 

	### ONEWAYSPL FUN ###
	# Function which calculates the Sound Pressure Level of a sound as it travels a 'one-way' route in a strait line of length R 

	onewaySPL<-function(sl,r,alpha,ts=0){
		sl - 20*log10(r)+ alpha*r +ts
	}

	#the SPL at which the sound hits the bat:

	preimpactSPL<-onewaySPL(SL_1m,preCDist,ALPHA)# SPL at which the sound impacts a bat 
	postimpactSPL<-onewaySPL(preimpactSPL,postCDist,ALPHA,batTS) # SPL at which the sound echoes off the bat

	hearTest = postimpactSPL >= hearTh
	return(hearTest)
}

rowcount = 0
for (i in 1:dim(cal)[2]){
# for each bat in the simulation
	rowcount = rowcount + dim(get(paste("calls.others",i-1,sep="")))[1]
	# increment rowcount with the number of rows contained in the agent's 
	# calls.others data frame
}	

echo = as.data.frame(matrix(NA, ncol = 5, nrow = rowcount))
# Empty data frame to store Echo data, of dimensions: 
# number of timesteps at which echoes are heard * 5
colnames(echo) = c("id.E","time.E","time.B","id.B","time.C")
# Columns: 
# - ID of the bat who emitted the call and obtained an echo in return
# - Time at which the call was emitted
# - Time at which the call's sound bounced off another bat
# - ID of the bat on which the sound bounced-off
# - Time at which the echo reached the call's emitter back.
meter = 0

for (i in 0:(dim(cal)[2]-1)){
# for each bat i in the simulation
	call.rec = get(paste("calls.others",i,sep=""))
	# get the corresponding data frame with heard calls from others
	firstcol = seq(1,dim(mov)[2],2)
	# store the index of the first column in mov data frame
	# with attributes of bat i
	
	for (j in 1:dim(call.rec)[1]){
	# for each call j heard by the bat i
		time.R = call.rec$Time_R[j]
		# time at which the call was heard
		# R stands for reception
		time.E = call.rec$Time_E[j]
		# time at which the call was emitted
		# E stands for emission
		id.E = call.rec$ID_E[j]
		# ID of the bat who emitted the call

		x.source = mov[time.R+1,firstcol[i+1]]
		# position of bat i on x-axis at the time it heard sound j
		y.source = mov[time.R+1,firstcol[i+1]+1]
		# position of bat i on y-axis at the time it heard sound j
		init.x = mov[time.E+1,firstcol[id.E+1]]
		# position of bat E on x-axis at the time call j
		# was heard by bat i.
		init.y = mov[time.E+1,firstcol[id.E+1]+1]
		# position of bat E on y-axis at the time call j
		# was heard by bat i.
		
		echoTravel = R.dist(0:(time.R - time.E + 5))
		preCD = (time.R - time.E)*VSOUND*TIMERESOLUTION
		postCD = echoTravel*VSOUND*TIMERESOLUTION
		echoPerception = Hear.echo(preCD, postCD, SOURCELEVEL, HEARINGTHRESHOLD)

		if (echoPerception){
			echo$time.C[meter+j] = time.R + echoTravel
			# time at which the collision C occurs, i.e. 
			# bat E hears back the echo bouncing-off from bat i
		}else{
			echo$time.C[meter+j] = NA
			# bat doesn't hear the echo because it is too faint
		}
				
		echo$id.B[meter+j] = i
		# ID of bat i, on which the call bounced-off 
		# & produced an echo
		echo$id.E[meter+j] = id.E
		# ID of bat E, who emitted the call 
		# & hears its echo back
		echo$time.B[meter+j] = time.R
		# time at which the call j bounced-off bat i
		echo$time.E[meter+j] = time.E
		# time at which the call was emitted

	}
	meter = meter + j		
}		


#---------- End of ECHOES CALCULATIONS ----------#

#---------- OVERLAP ANALYSIS ----------#
# This section aims to quantitatively analyse the interference between a bat's own echoes, 
# and the calls that are emitted by other bats and heard by the bat itself.
# To do so, we will look at:
# - mean number of calls/echo overlapping with an own-echo +/- SD
# - mean amount of time overlap each echo has with calls/echoes +/- SD
# - total %age of echoes that face overlapping with calls/echoes
# - the amount of time/ %age of the IPI that is free of any calls/echoes
# > P(overlap echo-echo) * P(overlap echo-call) --> how many echoes are distinctly heard at all
 
# These calculations will be made for both individual & group level,
# with overall estimations and description of the evolution of the estimators
# with IPI, Population size, Call Duration.

# Need to filter out sounds that reach the agent when it was calling, 
# rendering it impossible for the agent to hear it 
# NB: Filter applied to both calls & echoes 

plotRanges <- function(x, xlim = x, main = deparse(substitute(x)),
	col = "black", sep = 0.5, ...){
	height <- 1
	if (is(xlim, "Ranges")) xlim <- c(min(start(xlim)), max(end(xlim)))
	bins <- disjointBins(IRanges(start(x), end(x) + 1))
	plot.new()
	plot.window(xlim, c(0, max(bins)*(height + sep)))
	ybottom <- bins * (sep + height) - height
	rect(start(x)-0.5, ybottom, end(x)+0.5, ybottom + height, col = col, ...)
	title(main)
	axis(1)
}

npop = 0:(dim(cal)[2]-1)

ECoverlapMat = matrix(NA, nrow = dim(cal)[2], ncol = 7)
ECoverlapMat[,1] = npop
colnames(ECoverlapMat)=c("ID","meanNumOvPerEcho","sdNumOv","meanTimeOvPerEch","sdTimeOv","%echoOverlapped","%timeFreeIPI")

EEoverlapMat = matrix(NA, nrow = dim(cal)[2], ncol = 7)
EEoverlapMat[,1] = npop
colnames(EEoverlapMat)=colnames(ECoverlapMat)		

for (i in npop){
	echo_channel_nofilter = echo[which(echo$id.E == i),] 
	call_channel_nofilter = get(paste("calls.others",i,sep=""))
	filter_channel = get(paste("calls.self",i,sep=""))	
	c.ranges = IRanges()
	e.ranges = IRanges()
	f.ranges = IRanges()

	f.temp = rep(FALSE, SIMDURATION)
	filtertimes = filter_channel$Time_R+1
	f.temp[filtertimes] = TRUE
	f.ranges = append(f.ranges, IRanges(f.temp))
	f.ranges = shift(f.ranges, -1)


	for (j in 1:length(npop-1)){
		ident = npop[-(i+1)][j]
		
		c.temp = rep(FALSE, SIMDURATION)
		calltimes = call_channel_nofilter$Time_R[which(call_channel_nofilter$ID_E == ident)]+1
		checkcall = calltimes[which(calltimes<length(c.temp))]
		c.temp[checkcall] = TRUE
		c.ranges = append(c.ranges, IRanges(c.temp))

		e.temp = rep(FALSE, SIMDURATION)
		echotimes = echo_channel_nofilter$time.C[which(echo_channel_nofilter$id.B == ident)]+1
		checkecho = echotimes[which(echotimes<length(e.temp))]
		e.temp[checkecho] = TRUE
		e.ranges = append(e.ranges, IRanges(e.temp))
	}	

	c.ranges = shift(c.ranges, -1)
	e.ranges = shift(e.ranges, -1)

	c.ranges.f = IRanges()
	e.ranges.f = IRanges()	

	for (n in 1:length(c.ranges)) c.ranges.f = append(c.ranges.f,gaps(resize(f.ranges, f.ranges@width-1, fix="start"), c.ranges@start[n], c.ranges@start[n]+c.ranges@width[n]-1)) 
	for (n in 1:length(e.ranges)) e.ranges.f = append(e.ranges.f,gaps(resize(f.ranges, f.ranges@width-1, fix="start"), e.ranges@start[n], e.ranges@start[n]+e.ranges@width[n]-1)) 
	
	ECOverlaps = findOverlaps(e.ranges.f, c.ranges.f)
	EEOverlaps = findOverlaps(e.ranges.f, e.ranges.f)
	EEOverlaps = EEOverlaps[-which(EEOverlaps@from == EEOverlaps@to)]

	numECOverlaps = rep(NA, length(e.ranges.f)) 
	numEEOverlaps = rep(NA, length(e.ranges.f))	

	widthECOverlaps = rep(NA, length(e.ranges.f)) 
	widthEEOverlaps = rep(NA, length(e.ranges.f))

	all.ranges = IRanges(rep(TRUE, SIMDURATION))
	ipi.ranges = gaps(f.ranges, all.ranges@start, all.ranges@start+all.ranges@width-1)
	ipi.ranges = ipi.ranges + 1
	c.free.ipi = IRanges()
	e.free.ipi = IRanges()
	for (n in 1:length(ipi.ranges)){ 
		c.free.ipi = append(c.free.ipi, gaps(c.ranges, start=ipi.ranges[n]@start, end=ipi.ranges[n]@start+ipi.ranges[n]@width-1))
		e.free.ipi = append(e.free.ipi, gaps(e.ranges, start=ipi.ranges[n]@start, end=ipi.ranges[n]@start+ipi.ranges[n]@width-1))
	}
	c.free.ipi = c.free.ipi + 1
	e.free.ipi = e.free.ipi + 1

	for (k in 1:length(e.ranges.f)){
		indivEC = which(ECOverlaps@from == k)
		indivEE = which(EEOverlaps@from == k)

		if (length(indivEC)>0){

			Covrange = ranges(ECOverlaps[indivEC], e.ranges.f, c.ranges.f)
			eraseC = which(Covrange@start == Covrange@start+Covrange@width-1)
			if (length(eraseC)>0) Covrange = Covrange[-eraseC]

			if (length(Covrange)>0){
				widthovC = sum(reduce(Covrange)@width-1)
				widthECOverlaps[k] = widthovC/(e.ranges.f@width[k]-1)
				numECOverlaps[k] = length(Covrange)

			}else{
				widthECOverlaps[k] = 0
				numECOverlaps[k] = 0
			}
		
		}else{
			widthECOverlaps[k] = 0
			numECOverlaps[k] = 0
		}

		if (length(indivEE)>0){	
			
			Eovrange = ranges(EEOverlaps[indivEE], e.ranges.f, e.ranges.f)
			eraseE = which(Eovrange@start == Eovrange@start+Eovrange@width-1)
			if (length(eraseE)>0) Eovrange = Eovrange[-eraseE]

			if (length(Eovrange)>0){
				widthovE = sum(reduce(Eovrange)@width-1)
				widthEEOverlaps[k] = widthovE/(e.ranges.f@width[k]-1)
				numEEOverlaps[k] = length(Eovrange)

			}else{
				widthEEOverlaps[k] = 0
				numEEOverlaps[k] = 0
			}

		}else{
			widthEEOverlaps[k] = 0
			numEEOverlaps[k] = 0
		}
		
	}
	
	ECoverlapMat[i+1,2] = mean(numECOverlaps) # mean number of overlaps/echo
	ECoverlapMat[i+1,3] = sd(numECOverlaps) # standard deviation of number of overlaps/echo
	ECoverlapMat[i+1,4] = 100*mean(widthECOverlaps) # mean % time of overlap 
	ECoverlapMat[i+1,5] = 100*sd(widthECOverlaps) # standard deviation % time of overlap
	ECoverlapMat[i+1,6] = 100*sum(numECOverlaps>0)/length(e.ranges.f) # total % echoes (number) masked by others' calls 
	ECoverlapMat[i+1,7] = 100*sum(c.free.ipi@width-1)/sum(ipi.ranges@width-1) # total % of IPI time that is free of calls from others 
	
	EEoverlapMat[i+1,2] = mean(numEEOverlaps) # number of echo overlaps/echo 
	EEoverlapMat[i+1,3] = sd(numEEOverlaps) # standard deviation of number of echo overlaps/echo
	EEoverlapMat[i+1,4] = 100*mean(widthEEOverlaps) # mean % time/length of overlaps
	EEoverlapMat[i+1,5] = 100*sd(widthEEOverlaps) # standard deviation % time/length of overlaps
	EEoverlapMat[i+1,6] = 100*sum(numEEOverlaps>0)/length(e.ranges.f)# total % echoes (number) masked by echoes 
	EEoverlapMat[i+1,7] = 100*sum(e.free.ipi@width-1)/sum(ipi.ranges@width-1) # total % of IPI time that is free of echoes 

}

outputFile1 = paste(newDir,"/echo_call_index.csv", sep="")
outputFile2 = paste(newDir,"/echo_echo_index.csv", sep="")

write.csv(ECoverlapMat, file=outputFile1)
write.csv(EEoverlapMat, file=outputFile2)

}
}

print(Sys.time() - startTime)
