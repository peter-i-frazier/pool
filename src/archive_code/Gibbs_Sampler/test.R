rm(list=ls())
source('gibbs_util.R')
source('../Naive_Bayes/getFeatures.R')
## set parameters
C <- 38
Ntrain <- 50
Ntest <- 50
nVal <- 8
## create artificial theta for hit and nothit
alpha.hit <- runif(8) * 20
alpha.nothit <- abs(1 - alpha.hit + runif(8) * 25)
theta.hit <- t(rdirichlet(C, alpha.hit))
theta.nothit <- t(rdirichlet(C, alpha.nothit))

## create artifical training data and test data,
## 
train <- matrix(1, nrow=Ntrain, ncol=C)
for (r in 1:Ntrain) {
	pp <- runif(C)
	for (i in 1:8) {
		train[r,] <- train[r,] + ((pp-theta.hit[i,])>0)
		pp <- pp-theta.hit[i,]
	}
}

test <- matrix(1, nrow=Ntest, ncol=C)
for (r in 1:Ntest) {
	pp <- runif(C)
	for (i in 1:8) {
		test[r,] <- test[r,] + ((pp-theta.nothit[i,])>0)
		pp <- pp-theta.nothit[i,]
	}
}

# ##### test code 1
dataFile <- '../../data/binaryData_v2.csv'
classFile <- '../../data/Reduced_AA_Alphabet.csv'

data.org <- data.frame(read.csv(dataFile, header = T, as.is = T, sep = ","))
AAclass <- read.csv(classFile, header=T, as.is = T, sep=",")

nL <- 19
nR <- 19
trainData <- getFeatures(data.org,AAclass,nL,nR)
No.testcases <- 1000
Size.library <- 300
burnin.step <- 3000
record.step <- 500
X.train <- as.matrix(trainData[,c(1:(nL+nR))])
colnames(train) <- colnames(X.train)

train.data <- train
test.data <- rbind(train, test)
# test.data <- ceiling(matrix(nVal*runif(No.testcases * dim(hit.X)[2]), ncol=dim(hit.X)[2]))
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




##### test code 2
wholedata <- rbind(train, test)
prob.test2 <- rep(0,Ntrain+Ntest)
Y2 <- c(rep(1,Ntrain), rep(0,Ntest))
for (n in 1:500) {
	print (n)
	theta.comb2 <- sampleTheta(wholedata, Y2, nVal)
	theta.12 <- theta.comb2$theta_1
	theta.02 <- theta.comb2$theta_0
	prob.test2 <- prob.test2 + getProb(wholedata, theta.12, theta.02)
}
prob.test2 <- prob.test2/500


