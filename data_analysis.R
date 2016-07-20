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

###------------###
### PARAMETERS ###
###------------###

SIMDUR = 20
# NB: duration of the simulation must be known beforehand in order 
# for the script to work

###------------------------------###
### DATA IMPORTATION & RESHAPING ###
###------------------------------###

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
mov = as.data.frame(matrix(NA, ncol = length(M.fileNames), nrow = SIMDUR))
id_M = ceiling(1:length(M.fileNames)/2)
colnames(mov) = paste("Bat",id_M,c(".X",".Y"), sep = "")

for (i in id_M){
	mov[,1+2*(i-1)] = read.table(paste(M.path,M.fileNames[1+2*(i-1)], sep=""), header = F, sep = "\t", dec = ".")
	mov[,2+2*(i-1)] = read.table(paste(M.path,M.fileNames[2+2*(i-1)], sep=""), header = F, sep = "\t", dec = ".")
}

### Hearing ###

library(dplyr)

# Data format I'm trying to acheive:

# Sim.Step	Bat.R		Bat.E		t.E
# 0		0		0		0
# 0		1		1		0
# 1		1		0		0
# 1		1		1		0...

# Sim.Step is the time step,  
# Bat.R is the ID of the bat who heard a sound (R stands for receptor), 
# Bat.E is the ID of the bat who produced the sound (E stands for emitter),
# t.E is the time at which the sound was emitted

H.path = "Hearing/"
H.fileNames = dir(H.path, pattern =".txt") # alphabetical order, i.e.:
# file names are ordered in IDs ascending order, and then: c, i, t.

hea = as.data.frame(matrix(NA, nrow = 1, ncol = 4))
colnames(hea) = c("Sim.Step","Bat.R","Bat.E","t.E")

for (i in seq(1,length(H.fileNames),3)){ 
	t.E = read.table(paste(H.path,H.fileNames[i], sep=""), header = F, sep = "\t", dec = ".")
	Bat.E = read.table(paste(H.path,H.fileNames[i+1], sep=""), header = F, sep = "\t", dec = ".")
	Sim.Step = read.table(paste(H.path,H.fileNames[i+2], sep=""), header = F, sep = "\t", dec = ".")
	Bat.R = rep(as.numeric(substring(H.fileNames[i], 1,1)),dim(Sim.Step)[1])
	Bat.R = as.data.frame(Bat.R)
	temp = bind_cols("Sim.Step" = Sim.Step, "Bat.R" = Bat.R, "Bat.E" = Bat.E, "t.E" = t.E) 
	colnames(temp) = c("Sim.Step","Bat.R","Bat.E","t.E")
	hea = bind_rows(hea,temp)
}	 

hea = hea[-1,]

###-----------------------------###
### SONAR INTERFERENCE ANALYSIS ###
###-----------------------------###