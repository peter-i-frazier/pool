rm(list=ls())

source('gibbs_util.R')

data <- as.matrix(read.csv('feature.csv', header=F, as.is=T))
data <- data[-1,]
X.0 <- data[9:13,]
#leave one out cross validation

	# Initialization
	X <- data[1:8,]
	Z <- c(1:8)
	W <- X
	for (i in 1:500) {
		W <- rbind(W, ceiling(6 * runif(dim(W)[2])))
	}
	theta.1 <- matrix(1/6,nrow=6,ncol=dim(W)[2])
	theta.0 <- theta.1
	Y <- sampleY(theta.1, theta.0, W, Z)
	# Burn in step
	for (t in 1:5000) {
		theta.comb <- sampleTheta(W,Y)
		theta.1 <- theta.comb$t1
		theta.0 <- theta.comb$t0
		Y <- sampleY(theta.1, theta.0, W, Z)
		W <- sampleW(Y, theta.1, theta.0, X, Z)
		Z <- sampleZ(W, X)
	}
	# # Sample step
	# prob.1 <- rep(0,8)
	# prob.0 <- rep(0,5)
	# for (t in 1:500) {
	# 	theta.comb <- sampleTheta(W,Y)
	# 	theta.1 <- theta.comb$t1
	# 	theta.0 <- theta.comb$t0
	# 	Y <- sampleY(theta.1, theta.0, W, Z)
	# 	W <- sampleW(Y, theta.1, theta.0, X, Z)
	# 	Z <- sampleZ(W, X)
	# 	prob.1 <- prob.1 + getProb(X, theta.1, theta.0)
	# 	prob.0 <- prob.0 + getProb(X.0, theta.1, theta.0)
	# }
	# prob.1 <- prob.1/500
	# prob.0 <- prob.0/500
	
	# Sample step
	theta.1.ave <- matrix(0,nrow=dim(theta.1)[1], ncol=dim(theta.1)[2])
	theta.0.ave <- theta.1.ave
	for (t in 1:500) {
		theta.comb <- sampleTheta(W,Y)
		theta.1 <- theta.comb$t1
		theta.0 <- theta.comb$t0
		Y <- sampleY(theta.1, theta.0, W, Z)
		W <- sampleW(Y, theta.1, theta.0, X, Z)
		Z <- sampleZ(W, X)
		theta.1.ave <- theta.1.ave + theta.1
		theta.0.ave <- theta.0.ave + theta.0
	}
	theta.1.ave <- theta.1.ave/500
	theta.0.ave <- theta.0.ave/500
		prob.hit <- getProb(X, theta.1.ave, theta.0.ave)
		prob.nothit <- getProb(X.0, theta.1.ave, theta.0.ave)
	

	# theta.ave.11 <- theta.1
	# theta.ave.00 <- theta.0
	# for (t in 1:500) {
	# 	print (t+3500)
	# 	theta.comb <- sampleTheta(W,Y)
	# 	theta.1 <- theta.comb$t1
	# 	theta.0 <- theta.comb$t0
	# 	Y <- sampleY(theta.1, theta.0, W, Z)
	# 	W <- sampleW(Y, theta.1, theta.0, X, Z)
	# 	Z <- sampleZ(W, X)
	# 	theta.ave.11 <- theta.ave.11 + theta.1
	# 	theta.ave.00 <- theta.ave.00 + theta.0
	# }
	# theta.ave.11 <- theta.ave.11 /501
	# theta.ave.00 <- theta.ave.00 /501
