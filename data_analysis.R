### BAT IBM SENSORY JAMMING DATA ANALYSIS ###

### R Script
### Created on 22/06/2016 at 15:15
### Author: Claire

rm(list = ls())
# clear workspace
setwd("D:/Bat_Project/Res") 
# Change this to the directory where the specific results that you
# want to analyse are stored. Should be a directory in Res, and contain
# 3 directories: Calling, Moving, Hearing containing the corresponding 
# data.

### Import & reshape data

SIMDUR = 900
# duration of the simulation must be known beforehand 
# in order for the script to work

### 1) Calling

C.path = "Calling/"
C.fileNames = dir(C.path, pattern =".txt") # alphabetical order, i.e.:
# file names are ordered in IDs ascending order
cal = as.data.frame(matrix(NA, ncol = length(C.fileNames), nrow = SIMDUR))

for (i in 1:length(C.fileNames)){
	cal[i] = read.table(paste(C.path,C.fileNames[i], sep=""), header = F, sep = "\t", dec = ".")
	colnames(cal)[i] = paste("B_",i, sep="")
}

### 2) Moving

M.path = "Moving/"
M.fileNames = dir(M.path, pattern =".txt") # alphabetical order, i.e.:
# file names are ordered in IDs ascending order, x first and y second.
mov = as.data.frame(matrix(NA, ncol = length(M.fileNames), nrow = SIMDUR+1))
id_M = ceiling(1:length(M.fileNames)/2)

# WARNING: here we need SIMDUR+1 rows, but SIMDUR only for Calling... ??

for (i in id_M){
	for (j in seq(1, length(M.fileNames), 2)){
		for (k in seq(2, length(M.fileNames), 2)){
			mov[i] = read.table(paste(M.path,M.fileNames[i], sep=""), header = F, sep = "\t", dec = ".")
			colnames(mov)[j] = paste("B_",id_M[j],"_x", sep="")
			colnames(mov)[k] = paste("B_",id_M[k],"_y", sep="")
		}
	}
}

# WARNING: Pb in iterating id name in names of columns: 
# returns 2 all the time

### 3) Hearing
### Here, bats do not necessarily hear the same amount of sounds, so 
### there's gonna be some NA
### Data format I'm trying to acheive:

### Sim.Step	Bat1.Who	Bat1.When	Bat2.Who	Bat2.When
### 1			NA		NA		1		1
### 2			1		1		NA		NA	
### 3			2		1		2		1
### End		2		156		1		523

### Sim.Step is the time step, BatX.Who is the identity of the bat heard
### by the X-th individual, and BatX.When is the time of call emission.

H.path = "Hearing/"
H.fileNames = dir(H.path, pattern =".txt") # alphabetical order, i.e.:
# file names are ordered in IDs ascending order, and then: c, i, t.
hea = as.data.frame(matrix(NA, ncol = (2/3)*length(H.fileNames)+1, nrow = SIMDUR))
hea[,1] = 1:dim(hea)[1]
colnames(hea) = c("Tmstp",paste("Bat",ceiling(1:((2/3)*length(H.fileNames)/2)),c(".Who",".When"), sep = ""))

# WARNING: colnames does not work yet!

for (i in 1:length(H.fileNames)){
	hea = merge(hea, read.table(paste(C.path,C.fileNames[i], sep=""), header = F, sep = "\t", dec = ".")
}










