###################################
#####                         #####
##### BAT IBM SENSORY JAMMING ##### 
#####                         #####
###################################

### Created on 22/06/2016 at 15:15
### Author: Claire

rm(list = ls())
# clear workspace
setwd("D:/Bat_Project/Res")
# Change this to the directory where the specific results that you want
# to analyse are stored. Should be a directory in Res, and contain 3
# directories: Calling, Moving, Hearing containing the corresponding data

SIMDURATION <- 20
# as.numeric(readline("Duration of simulation: ")) # iterations
TIMERESOLUTION <- 0.002
# as.numeric(readline("Time resolution of the simulation: ")) # seconds per iteration
FLIGHTSPEED <- 5.5 
# as.numeric(readline("Velocity of flight for all bats: ")) # m/s
VSOUND <- 340.29 # m/s
TETA = 0 # bat movement angle, in radians

### Calling ###

C.path = "Calling/"
C.fileNames = dir(C.path, pattern =".txt") # alphabetical order, i.e.:
# file names are ordered in IDs ascending order
cal = as.data.frame(matrix(NA, ncol = length(C.fileNames), nrow = SIMDURATION))

for (i in 1:length(C.fileNames)){
	cal[,i] = read.table(paste(C.path,C.fileNames[i], sep=""), header = F, sep = "\t", dec = ".")
	colnames(cal)[i] = paste("B_",i, sep="")
}

### Moving ###

M.path = "Moving/"
M.fileNames = dir(M.path, pattern =".txt") # alphabetical order, i.e.:
# file names are ordered in IDs ascending order, x first and y second.
mov = as.data.frame(matrix(NA, ncol = length(M.fileNames), nrow = SIMDURATION))
id_M = ceiling(1:length(M.fileNames)/2)
colnames(mov) = paste("Bat",id_M,c(".X",".Y"), sep = "")

for (i in id_M){
	mov[,1+2*(i-1)] = read.table(paste(M.path,M.fileNames[1+2*(i-1)], sep=""), header = F, sep = "\t", dec = ".")
	mov[,2+2*(i-1)] = read.table(paste(M.path,M.fileNames[2+2*(i-1)], sep=""), header = F, sep = "\t", dec = ".")
}

# From this, we extract the inter-individual distances...
batPop = (1:(dim(mov)[2]/2))
Column = seq(1,dim(mov)[2], 2)
x = as.vector(rep(NA,dim(mov)[2]/2))
y = as.vector(rep(NA,dim(mov)[2]/2))

for (i in batPop){
	x[i] = mov[1,Column[i]]
	y[i] = mov[1,Column[i]+1]
}

IID = dist(cbind(x,y)) # Matrix of IID

### Hearing ###

# Data format I'm trying to acheive:

# Sim.Step	Bat.R		Bat.E		t_E
# 0		0		0		0
# 0		1		1		0
# 1		1		0		0
# 1		1		1		0...

# Sim.Step is the time step,  
# Bat.R is the ID of the bat who heard a sound (R stands for receptor, 
# Bat.E is the ID of the bat who produced the sound (E stands for emitter),
# t_E is the time at which the sound was emitted

H.path = "Hearing/"
H.fileNames = dir(H.path, pattern =".txt") # alphabetical order, i.e.:
# file names are ordered in IDs ascending order, and then: c, i, t.

for (i in seq(1,length(H.fileNames),3)){
   	
	ownID = (ceiling(1:length(H.fileNames)/3)-1)[i]
	tcall = read.table(paste(H.path,H.fileNames[i], sep=""), header = F, sep = "\t", dec = ".")
	idbat = read.table(paste(H.path,H.fileNames[i+1], sep=""), header = F, sep = "\t", dec = ".")
	tmstp = read.table(paste(H.path,H.fileNames[i+2], sep=""), header = F, sep = "\t", dec = ".")

	temp = as.data.frame(matrix(NA, ncol = 3, nrow = dim(tcall)[1]))
	colnames(temp) = c("Time_R","ID_E","Time_E")
	temp$Time_R = tmstp[,1]
	temp$ID_E = idbat[,1]
	temp$Time_E = tcall[,1]
	own_temp = temp[which(temp$ID_E == ownID),]
	other_temp = temp[-which(temp$ID_E == ownID),]
	D_E_R = (other_temp$Time_R - other_temp$Time_E)*TIMERESOLUTION*VSOUND
	other_temp_bis = cbind(other_temp, D_E_R)
	assign(paste("calls.self",ownID,sep=""),own_temp)
	assign(paste("calls.others",ownID,sep=""),other_temp)
}

# Echoes from Time_R and Time_E
# Calculate IID
# Calculate distance between call source of any other bat and position where focal bat registered the sound
# Calculate distance travelled by bats in the mean time

R.dist = function(expl.t){
	new.x = init.x + FLIGHTSPEED * expl.t * cos(TETA)
	new.y = init.y + FLIGHTSPEED * expl.t * sin(TETA)
	R.t.sq = (new.x - X.SOURCE) ** 2 + (new.y - Y.SOURCE) ** 2
	r.t.sq = (VSOUND * expl.t) ** 2
		
	collision = which(R.t.sq <= r.t.sq)
	impact.t = collision[length(collision)]		
	return(impact.t)
}

echo = as.data.frame(matrix(NA, ncol = 5, nrow = dim(call.rec)[1] * dim(cal)[2]))
colnames(echo) = c("id.E","time.E","time.C","id.C","time.D")
meter = 0

for (i in 0:(dim(cal)[2]-1)){
	call.rec = get(paste("calls.others",i,sep="")) 

	for (j in 1:dim(call.rec)[1]){
		firstcol = seq(1,dim(cal)[2]*2,2)
		time.R = call.rec$Time_R[j]
		time.E = call.rec$Time_E[j]
		id.em = call.rec$ID_E[j]

		X.SOURCE = mov[firstcol[i+1],j]
		Y.SOURCE = mov[firstcol[i+2],j]
		init.x = mov[firstcol[id.em],j]
		init.y = mov[firstcol[id.em+1],j]
		
		echo$time.D[meter+j] = R.dist(0:(time.R - time.E))
		echo$id.C[meter+j] = i
		echo$id.E[meter+j] = id.em
		echo$time.C[meter+j] = time.R
		echo$time.E[meter+j] = time.E

	}
	meter = meter + j		
}		




