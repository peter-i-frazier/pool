rm(list=ls())
library(doMC)
library(foreach)
source('sparse_prior_util.R')
#=================================================================================
#Specify Paths and working directory
dataFile <- '../../data/binaryData.csv'
classFile <- '../../data/Reduced_AA_Alphabet.csv'

#get data 
data.org <- data.frame(read.csv(dataFile, header = T, as.is = T, sep = ","))
AAclass <- read.csv(classFile, header=T, as.is = T, sep=",")

#=================================================================================
#Set parameters
nL <- 19
nR <- 19
trainData <- getFeatures(data.org,AAclass,nL,nR)
nAA <- max(AAclass[1,])
burnin.step <- 5000
record.step <- 1000
#=================================================================================
## For AcpH
outcome <- 'AcpH'
X <- as.matrix(trainData[,c(1:(nL+nR))])
Y <- trainData[,outcome]

## cross validation
registerDoMC()
prob <- foreach (n=1:length(Y), .combine=rbind) %dopar% {
	sparsePrior(X[-n,], Y[-n], X[n,], nAA, burnin.step, record.step)
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