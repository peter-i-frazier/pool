rm(list=ls())
library(doMC)
library(foreach)
#=================================================================================
#Specify Paths and working directory
dataPath <- 'D:/Study/Summer2013/Peptide/data.original'
dataFile <- paste(dataPath,  '/newData.csv', sep = "")
classFile <- paste(dataPath, '/Reduced_AA_Alphabet.csv', sep = "")
srcPath <- 'D:/Study/Summer2013/Peptide/src'
wrkDir <- 'D:/Study/Summer2013/Peptide/wrkDir'
setwd(wrkDir)

source(paste(srcPath, '/Sparse_prior/sparse_prior_util.R', sep = ""))

#get data 
data_org <- data.frame(read.csv(dataFile, header = T, as.is = T, sep = ","))
AAclass <- read.csv(classFile, header=T, as.is = T, sep=",")
#=================================================================================
#Parallel Computation Packages
require(multicore)
require(snow)
require(doParallel)
require(MCMCpack)
require(Rlab)
#=================================================================================
#Parallel Computation Setting
cl <- makeCluster(4)
registerDoParallel(cl)
#=================================================================================
#Set parameters
nL <- 19
nR <- 19
trainData <- getFeatures(data_org,AAclass,nL,nR)
nAA <- max(AAclass[1,])
burnin.step <- 5000
record.step <- 10000
#=================================================================================
## For AcpH
outcome <- 'AcpH'
X <- as.matrix(trainData[,c(1:(nL+nR))])
Y <- trainData[,outcome]
tX <- X[Y==0,][-c(1:6),]

Sample <- sample.int(dim(tX)[1], size = dim(X[Y==1,])[1]-6)
trainX <- rbind(X[c(1:16),], X[Y==1,][-c(1:10),], tX[Sample,])
trainY <- c(rep(1,10), rep(0,6), rep(1,dim(X[Y==1,])[1]-10),rep(0,dim(X[Y==1,])[1]-6))



## cross validation
registerDoMC()
prob_1 <- foreach (n = 1:16, .init = c(), .combine = rbind, .packages = c('MCMCpack','Rlab'), .inorder = TRUE, .errorhandling = 'pass') %dopar% {
	sparsePrior(trainX[-n,], trainY[-n], trainX[n,], nAA, burnin.step, record.step)
}

prob <- c()
for ( n in 1:16 ) {
	prob <- rbind(prob, sparsePrior(X[-n,], Y[-n], X[n,], nAA, burnin.step, record.step))
}
	
# test <- X[1,]
# X <- X[-1,]
# Y <- Y[-1]
# # some constant
# nF <- dim(X)[2]
# X.1 <- X[Y==1,]
# X.0 <- X[Y==0,]
# factor <- length(Y==1)/length(Y==0)
# # Initialization
# ztable <- Ztable(nAA)
# P.1 <- rep(0.5, nF)
# P.0 <- P.1
# Z.1 <- c()
# for (i in 1:nF) {
# 	Z.1 <- cbind(Z.1,rbern(nAA, P.1[i]))
# }
# Z.0 <- Z.1