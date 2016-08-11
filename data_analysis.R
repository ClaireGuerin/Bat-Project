###############################################
###                         				###
### 		BAT IBM SENSORY JAMMING 		### 
###                         				###
###############################################

### Created on 22/06/2016 at 15:15
### Author: Claire

rm(list = ls())
# clear current workspace
setwd("D:/Bat_Project/Res")
# Change this to the directory where the specific results that you want
# to analyse are stored. Should be a folder containing 3 sub-folders:
# Calling, Moving, Hearing, with the corresponding data

#---------- PARAMETERS OF THE SIMULATION ----------#
# This section should be modified for each data, extracted from
# a new simulation, that is to be analysed. If the parameters below 
# are not entered accordingly to the simulation parameters
# corresponding to the data, the script will not work.

SIMDURATION <- 30 # duration of the simulation iterations
TIMERESOLUTION <- 0.002 # time resolution, in seconds per iteration
FLIGHTSPEED <- 5.5 # speed of flight of the bats, in m/s
VSOUND <- 340.29 # m/s
TETA <- 0 # bat movement angle, in radians
SOURCELEVEL <- 120 # source level of the bat's call in dB SPL at 10 cm
HEARINGTHRESHOLD <- 0 # hearing threshold of the bat in dB SPL
ALPHA <- -1.7 # sound absorbtion rate at particular frequency in db/m  

#---------- End of SIMULATION PARAMETERS ----------#


#---------- DATA IMPORT & RESHAPE ----------#
# This section imports the data from the simulation 
# ran under python with the sensory_jamming.py program.
# It also reshapes the data into a R-friendly data format.

### CALLING ###

C.path = "Calling/" # set path to the Calling folder.
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

M.path = "Moving/" # set path to the Moving folder.
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

H.path = "Hearing/" # set path to the Hearing folder.
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

	new.x = init.x + FLIGHTSPEED * expl.t * cos(TETA)
	# new position of the moving bat on the x-axis, after delta t = expl.t time
	new.y = init.y + FLIGHTSPEED * expl.t * sin(TETA)
	# new position of the moving bat on the y-axis, after delta t = expl.t time
	R.t.sq = (new.x - x.source) ** 2 + (new.y - y.source) ** 2
	# R.t.sq stands for R(t)². Distance between the moving bat, which
	# originally emitted the sound, and the position where the echo 
	# bounced-off, otherwise called "source" as in source of the echo.
	# R(t)² slowly decreases, then increases with time.
	r.t.sq = (VSOUND * expl.t) ** 2
	# r.t.sq stands for r(t)². Radius of the sound as it propagates,
	# i.e. distance between the source of the echo and its current
	# position in space. r(t)² increases with time.
		
	collision = which(R.t.sq <= r.t.sq)
	# a bat-sound collision occurs as soon as R(t)² <= r(t)²
	impact.t = collision[1]
	# Exact time of impact
	return(impact.t)
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
	SL_1m = SL - abs(20*log10(0.1/1)) # Sound Level calculated at 1 m 

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
		init.x = mov[time.R+1,firstcol[id.E+1]]
		# position of bat E on x-axis at the time call j
		# was heard by bat i.
		init.y = mov[time.R+1,firstcol[id.E+1]+1]
		# position of bat E on y-axis at the time call j
		# was heard by bat i.
		
		echoTravel = R.dist(0:(time.R - time.E))
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
# - %age of echoes that face overlapping with calls/echoes
# - the amount of time/ %age of the IPI that is free of any calls/echoes
# > P(overlap echo-echo) * P(overlap echo-call) --> how many echoes are distinctly heard at all
 
# These calculations will be made for both individual & group level,
# with overall estimations and description of the evolution of the estimators
# with IPI, Population size, Call Duration.

# First, filter out sounds that reach the agent when it was calling, 
# rendering it impossible for the agent to hear it 
# NB: Filter applied to both calls & echoes 

for (i in 0:(dim(cal)[2]-1)){
	e.self = echo[which(echo$id.E == i),] 
	c.self = get(paste("calls.self",i,sep=""))
	c.others = get(paste("calls.others",i,sep=""))	
	
	for (j in 1:dim(c.self)[1]){
		filter1 = which(c.others$Time_R == c.self$Time_R[j])
		if (length(filter1) > 0){
			c.others = c.others[-filter1,]
		}
		
		filter2 = which(e.self$time.C == c.self$Time_R[j])
		if (length(filter2) > 0){
			e.self = e.self[-filter2,]		}		
	
	}
	name1 = paste("calls.filter",i,sep="")
	assign(name1, c.others[order(c.others$ID_E),])
	name2 = paste("echo.filter",i,sep="")
	assign(name2, e.self)
}

library(IRanges)

npop = 0:(dim(cal)[2]-1)

ECoverlapMat = matrix(NA, nrow = dim(cal)[2], ncol = 5)
ECoverlapMat[,1] = npop

EEoverlapMat = matrix(NA, nrow = dim(cal)[2], ncol = 5)
EEoverlapMat[,1] = npop

for (i in npop){
	echo_channel = get(paste("echo.filter",i,sep=""))
	call_channel = get(paste("calls.filter",i,sep=""))	
	c.ranges = IRanges()
	e.ranges = IRanges()

	for (j in 1:3){
		ident = npop[-(i+1)][j]
		
		c.temp = rep(FALSE, SIMDURATION)
		calltimes = call_channel$Time_R[which(call_channel$ID_E == ident)]+1
		c.temp[calltimes] = TRUE
		c.ranges = append(c.ranges, IRanges(c.temp))

		e.temp = rep(FALSE, SIMDURATION)
		echotimes = echo_channel$time.C[which(echo_channel$id.B == ident)]+1
		e.temp[echotimes] = TRUE
		e.ranges = append(e.ranges, IRanges(e.temp))
	}	

	c.ranges = shift(c.ranges, -1)
	e.ranges = shift(e.ranges, -1)
	
	ECOverlaps = findOverlaps(e.ranges, c.ranges)
	ECOverlapWidth = ranges(ECOverlaps, e.ranges, c.ranges)

	EEOverlaps = findOverlaps(e.ranges, e.ranges)
	EEOverlapWidth = ranges(EEOverlaps, e.ranges, e.ranges)

	numECOverlaps = rep(NA, length(e.ranges)) 
	numEEOverlaps = rep(NA, length(e.ranges))	

	for (k in 1:length(e.ranges)){
		numECOverlaps[k] = length(which(ECOverlaps@from == k))
		numEEOverlaps[k] = length(which(EEOverlaps@from == k))
	}
	
	ECoverlapMat[i+1,2] = mean(numECOverlaps) # mean number of overlaps/echo
	ECoverlapMat[i+1,3] = sd(numECOverlaps) # standard deviation of number of overlaps/echo
	ECoverlapMat[i+1,4] = mean(width(ECOverlapWidth)*TIMERESOLUTION) # mean time of overlap 
	ECoverlapMat[i+1,5] = sd(width(ECOverlapWidth)*TIMERESOLUTION) # standard deviation time of overlap 
	
	EEoverlapMat[i+1,2] = mean(numEEOverlaps) # number of echo overlaps/echo 
	EEoverlapMat[i+1,3] = sd(numEEOverlaps) # standard deviation of number of echo overlaps/echo
	EEoverlapMat[i+1,4] = mean(width(EEOverlapWidth)*TIMERESOLUTION) # mean time/length of overlaps
	EEoverlapMat[i+1,5] = sd(width(EEOverlapWidth)*TIMERESOLUTION) # standard deviation time/length of overlaps
	
}

EEoverlapMat
subsetByOverlaps(cr0,er0)
mergeOv=mergeByOverlaps(cr0,er0)
intersect(mergeOv[3,1],mergeOv[3,2])

Overlaps@queryHits # doesn't work for me O.o
