rm(list=ls())

source('gibbs_util.R')

data <- as.matrix(read.csv('feature.csv', header=F, as.is=T))
data <- data[-1,]
X.0 <- data[9:13,]
#leave one out cross validation
for (N in 1:8) {
	print (N)
	# Initialization
	X <- data[1:8,]
	X <- X[-N,]
	Z <- c(1:7)
	W <- X
	for (i in 1:500) {
		W <- rbind(W, ceiling(6 * runif(dim(W)[2])))
	}
	theta.1 <- matrix(1/6,nrow=6,ncol=dim(W)[2])
	theta.0 <- theta.1
	Y <- sampleY(theta.1, theta.0, W, Z)
	output <- list(Y=Y, W=W, X=X, Z=Z, t1=theta.1, t0=theta.0)
	# Burn in step
	output <- gibbsIter(output$Y, output$W, output$X, output$Z, 5000)
	# Sample step
	prob.1 <- rep(0,7)
	prob.0 <- rep(0,5)
	for (t in 1:500) {
		output <- gibbsIter(output$Y, output$W, output$X, output$Z)
		prob.1 <- prob.1 + getProb(X, output$t1, output$t0)
		prob.0 <- prob.0 + getProb(X.0, output$t1, output$t0)
	}
	prob.1 <- prob.1/500
	prob.0 <- prob.0/500

	# Combine results
	if (N ==1) {
		prob.hit <- prob.1
		prob.nothit <- prob.0
	}
	else {
		prob.hit <- cbind(prob.hit, prob.1)
		prob.nothit <- cbind(prob.nothit, prob.0)
	}
}
write.csv(prob.hit,'result/prob_hit_enzyme2.csv')
write.csv(prob.nothit,'result/prob_nothit_enzyme2.csv')
