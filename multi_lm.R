###############################################
###                         				###
### 		BAT IBM SENSORY JAMMING 		###
###		INFERENTIAL STATISTICS			### 
###                         				###
###############################################

### Created on 12/08/2016 at 11:33
### Author: Claire Guérin

rm(list = ls())
# clear current workspace

all_ipi = seq(0,1,0.01)
all_nedge = seq(2,10,1)
all_iidonaxe = seq(1,10,1)

all_comb = expand.grid(all_ipi,all_nedge,all_iidonaxe)
colnames(all_comb) = c("all_ipi","all_nedge","all_iidonaxe")
maxpop = max(all_iidonaxe)^2

indiv = as.data.frame(matrix(NA, nrow = dim(all_comb)[1], ncol = maxpop))
colnames(indiv) = as.character(1:maxpop-1)
data = cbind(all_comb, indiv)

data.expanded = data[rep(row.names(data), 6),]
rownames(data.expanded) = as.character(1:dim(data.expanded)[1])

###----------IMPORT INDICES----------###

setwd("D:/Bat_Project/Res/Overlap_indices")
fileNames = list.files() # test: list.files("D:/Bat_Project/Res/Calling")

# test: x="echo_call_ipi_nedge_iid.csv"

combInd = ceiling(1:dim(data.expanded)[1]/6)[1:18]
indexInd = rep(1:6, dim(data)[1])

for (i in 1:dim(data.expanded)[1]){
	combNum = combInd[i]
	indexNum = indexInd[i]
	comb = paste(data[combNum,1],"_",data[combNum,2],"_",data[combNum,3], sep="") #comb="_"
	pos = grep(comb, fileNames) # should give 2 files
	ec.file = grep("echo_call", fileNames[pos]) 
	ee.file = grep("echo_echo", fileNames[pos])
	
	ec.index = read.csv(fileNames[ec.file])
	ee.index = read.csv(fileNames[ee.file])

	index_indiv = ec.index[,indexNum+1]
	data.expanded[i,1:length(index_indiv)+3] = index_indiv
}

data.split = split(data.expanded,rep(1:dim(data)[1],each=6))
listDf = names(data.split)
sapply(

for (i in 1:length(data.split)){
	
	subDf = get(as.character(i),rdmat2)#data.split)
	subDf.T = t(subDf)
	x$subDf.T)
}

listDf = as.character(1:dim(data)[1])


# rdmat2 = split(rdmat, rep(1:2,2))
# test: rdmat=as.data.frame(matrix(NA, 4, 7))
# test: colnames(rdmat)=c("ID","meanNumOv","sdNumOv","meanTimeOv","sdTimeOv","numMask","timeFreeIPI")
# test: for (n in 1:dim(rdmat)[1]) rdmat[n,] = runif(7, 0.0, 1.0)
# test: h=colMeans(rdmat)[2:7]
# test: rdmat[1,2:7]=h
###----------End of INDICES IMPORT----------###

# test: set.seed(69)
# test: mockIndex = runif(dim(all_comb)[1], 0.0, 1.0)
# test: length(mockIndex)
# test: all_comb2 = cbind(all_comb, mockIndex)


###----------MULTIPLE LINEAR REGRESSIONS----------###

for (i in 1:6)
model1 = lm(mockIndex~all_ipi*all_nedge*all_iidonaxe, data = all_comb2)
summary(model1)
x=resid(model1)
hist(x) 
# no shapiro.test(x) possible because when the sample size is large, 
# even the smallest deviation from normality will lead to a rejected null.

library(MASS)
step=stepAIC(model1, direction="both")
step$anova # gives you the final model

library(scatterplot3d)
scatterplot3d(y = all_comb2$mockIndex, x = all_comb2$all_ipi, z = all_comb2$all_nedge)