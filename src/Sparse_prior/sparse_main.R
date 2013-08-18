rm(list=ls())
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
burnin.step <- 4000
record.step <- 500
#=================================================================================
## For AcpH
outcome <- 'AcpH'
X <- as.matrix(trainData[,c(1:(nL+nR))])
Y <- trainData[,outcome]

# ## cross validation
# prob <- c()
# for (n in 1:length(Y)) {
# 	print (n)
# 	X.train <- X[-n,]
# 	Y.train <- Y[-n]
# 	test <- X[n,]
# 	prob <- c(prob, mean(sparsePrior(X.train, Y.train, test, nAA, burnin.step, record.step)))
# }

test <- X
# some constant
nF <- dim(X)[2]
X.1 <- X[Y==1,]
X.0 <- X[Y==0,]
factor <- length(Y==1)/length(Y==0)
# Initialization
ztable <- Ztable(nAA)
P.1 <- rep(0.5, nF)
P.0 <- P.1
Z.1 <- c()
for (i in 1:nF) {
	Z.1 <- cbind(Z.1,rbern(nAA, P.1[i]))
}
Z.0 <- Z.1