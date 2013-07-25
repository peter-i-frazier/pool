rm(list=ls())
source('gibbs_util.R')
source('../Naive_Bayes/getFeatures.R')
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
nVal <- max(AAclass[1,])
No.testcases <- 1000
Size.library <- 500
burnin.step <- 4000
record.step <- 500
#=================================================================================
## For AcpH
outcome <- 'AcpH'
X.train <- as.matrix(trainData[,c(1:(nL+nR))])
Y.train <- trainData[,outcome]
hit.X <- X.train[Y.train==1,]
nothit.X <- X.train[Y.train==0,]

# # leave one out cross validation
# prob.hit <- c()
# N <- dim(hit.X)[1]
# for (n in 1:N) {
# 	print (n)
# 	# Initialization
# 	train <- hit.X
# 	test <- train[n,]
# 	train <- train[-n,]
# 	prob.1 <- gibbsSampler(train, test, AAclass, Size.library, burnin.step, record.step)
# 	prob.hit <- c(prob.hit, prob.1)
# }
# take in all the data and predict prob not hit
# prob.0 <- gibbsSampler(hit.X, nothit.X, nVal, Size.library, burnin.step, record.step)

# write.csv(prob.hit,'../result/prob_hit_AcpH.csv')
# write.csv(prob.0,'../result/prob_nothit_AcpH.csv')

# #=================================================================================
# ## For sfp
# outcome <- 'sfp'
# X.train <- as.matrix(trainData[,c(1:(nL+nR))])
# Y.train <- trainData[,outcome]
# hit.X <- X.train[Y.train==1,]

# # leave one out cross validation
# prob.hit <- c()
# N <- dim(hit.X)[1]
# for (n in 1:N) {
# 	print (n)
# 	# Initialization
# 	train <- hit.X
# 	test <- train[n,]
# 	train <- train[-n,]
# 	prob.1 <- gibbsSampler(train, test, AAclass, Size.library, burnin.step, record.step)
# 	prob.hit <- c(prob.hit, prob.1)
# }
# # Test random peptides' prob being hit
# test <- ceiling(matrix(max(AAclass[1,])*runif(No.testcases * dim(hit.X)[2]), ncol=dim(hit.X)[2]))
# prob.0 <- gibbsSampler(hit.X, test, AAclass, Size.library, burnin.step, record.step)

# write.csv(prob.hit,'../result/prob_hit_sfp.csv')
# write.csv(prob.0,'../result/prob_random_nothit_sfp.csv')



# test part
# some constant
	train.data <- hit.X
	# test.data <- ceiling(matrix(nVal*runif(No.testcases * dim(hit.X)[2]), ncol=dim(hit.X)[2]))
	test.data <- nothit.X
	M <- Size.library

	N <- dim(train.data)[1]
	C <- dim(train.data)[2]
	# Initialization
	Z <- c(1:N)
	W <- train.data
	for (i in 1:(M-N)) {
		W <- rbind(W, ceiling(nVal * runif(C)))
	}
	theta.1 <- matrix(1/nVal,nrow=nVal,ncol=C)
	theta.0 <- theta.1
	Y <- sampleY(theta.1, theta.0, W, Z)
	# Burn in step
	for (t in 1:burnin.step) {
		print (t)
		theta.comb <- sampleTheta(W, Y, nVal)
		theta.1 <- theta.comb$theta_1
		theta.0 <- theta.comb$theta_0
		Y <- sampleY(theta.1, theta.0, W, Z)
		W <- sampleW(Y, theta.1, theta.0, train.data, Z)
		Z <- sampleZ(W, train.data)
	}
	# Record step
	if (is.vector(test.data)) {
		prob <- 0
	} else {
		prob <- rep(0, dim(test.data)[1])
	}
	for (t in 1:record.step) {
		print (t)
		theta.comb <- sampleTheta(W, Y, nVal)
		theta.1 <- theta.comb$theta_1
		theta.0 <- theta.comb$theta_0
		Y <- sampleY(theta.1, theta.0, W, Z)
		W <- sampleW(Y, theta.1, theta.0, train.data, Z)
		Z <- sampleZ(W, train.data)
		prob <- prob + getProb(test.data, theta.1, theta.0)
	}
	prob <- prob/record.step

