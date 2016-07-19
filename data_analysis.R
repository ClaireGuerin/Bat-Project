###################################
#####                         #####
##### BAT IBM SENSORY JAMMING ##### 
#####                         #####
###################################

### Created on 22/06/2016 at 15:15
### Author: Claire

rm(list = ls())
# clear workspace
setwd("C:/Users/dlewanzik/Documents/Bat-Project/Res") 
# Change this to the directory where the specific results that you want
# to analyse are stored. Should be a directory in Res, and contain 3
# directories: Calling, Moving, Hearing containing the corresponding data

###------------------------------###
### DATA IMPORTATION & RESHAPING ###
###------------------------------###

SIMDUR = 900
# NB: duration of the simulation must be known beforehand in order 
# for the script to work

### Calling ###

C.path = "Calling/"
C.fileNames = dir(C.path, pattern =".txt") # alphabetical order, i.e.:
# file names are ordered in IDs ascending order
cal = as.data.frame(matrix(NA, ncol = length(C.fileNames), nrow = SIMDUR))

for (i in 1:length(C.fileNames)){
	cal[,i] = read.table(paste(C.path,C.fileNames[i], sep=""), header = F, sep = "\t", dec = ".")
	colnames(cal)[i] = paste("B_",i, sep="")
}

### Moving ###

M.path = "Moving/"
M.fileNames = dir(M.path, pattern =".txt") # alphabetical order, i.e.:
# file names are ordered in IDs ascending order, x first and y second.
mov = as.data.frame(matrix(NA, ncol = length(M.fileNames), nrow = SIMDUR+1))
id_M = ceiling(1:length(M.fileNames)/2)
colnames(mov) = paste("Bat",id_M,c(".X",".Y"), sep = "")

for (i in id_M){
	mov[,1+2*(i-1)] = read.table(paste(M.path,M.fileNames[1+2*(i-1)], sep=""), header = F, sep = "\t", dec = ".")
	mov[,2+2*(i-1)] = read.table(paste(M.path,M.fileNames[2+2*(i-1)], sep=""), header = F, sep = "\t", dec = ".")
}

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
	tcall = read.table(paste(H.path,H.fileNames[i], sep=""), header = F, sep = "\t", dec = ".")
	idbat = read.table(paste(H.path,H.fileNames[i+1], sep=""), header = F, sep = "\t", dec = ".")
	tmstp = read.table(paste(H.path,H.fileNames[i+2], sep=""), header = F, sep = "\t", dec = ".")

	for (j in 1:dim(tmstp)[1]){
		hea[which(hea$Tmstp == tmstp[j,1]),2+2*((i-1)/3)] = idbat[j,1]
		hea[which(hea$Tmstp == tmstp[j,1]),2+2*((i-1)/3)+1] = tcall[j,1]
		# hea[,2+2*((i-1)/3)] = as.factor(hea[,2+2*((i-1)/3)])	 


hea = as.data.frame(matrix(NA, ncol = (2/3)*length(H.fileNames)+1, nrow = SIMDUR))
hea[,1] = 0:(dim(hea)[1]-1)
colnames(hea) = c("Tmstp",paste("Bat",ceiling(1:((2/3)*length(H.fileNames))/2),c(".Who",".When"), sep = ""))

for (i in seq(1,length(H.fileNames),3)){ 
	tcall = read.table(paste(H.path,H.fileNames[i], sep=""), header = F, sep = "\t", dec = ".")
	idbat = read.table(paste(H.path,H.fileNames[i+1], sep=""), header = F, sep = "\t", dec = ".")
	tmstp = read.table(paste(H.path,H.fileNames[i+2], sep=""), header = F, sep = "\t", dec = ".")

	for (j in 1:dim(tmstp)[1]){
		hea[which(hea$Tmstp == tmstp[j,1]),2+2*((i-1)/3)] = idbat[j,1]
		hea[which(hea$Tmstp == tmstp[j,1]),2+2*((i-1)/3)+1] = tcall[j,1]
		# hea[,2+2*((i-1)/3)] = as.factor(hea[,2+2*((i-1)/3)])
	}
}

