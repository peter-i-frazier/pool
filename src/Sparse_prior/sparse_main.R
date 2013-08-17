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
burnin.step <- 50
record.step <- 50
#=================================================================================
## For AcpH
outcome <- 'AcpH'
X <- as.matrix(trainData[,c(1:(nL+nR))])
Y <- trainData[,outcome]

## cross validation
prob <- c()
for (n in 1:length(Y)) {
	print (n)
	X.train <- X[-n,]
	Y.train <- Y[-n]
	test <- X[n,]
	prob <- c(prob, mean(sparsePrior(X.train, Y.train, test, nAA, burnin.step, record.step)))
}

