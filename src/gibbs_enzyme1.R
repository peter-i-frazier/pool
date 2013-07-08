rm(list=ls())

source('gibbs_util.R')

data <- as.matrix(read.csv('feature.csv', header=F, as.is=T))
data <- data[-1,]
### leave one out cross validation
for (N in 1:13) {
	print (N)
	# Initialization
	X <- data
	X <- X[-N,]
	Z <- c(1:12)
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
	prob.1 <- rep(0,12)
	for (t in 1:500) {
		output <- gibbsIter(output$Y, output$W, output$X, output$Z)
		prob.1 <- prob.1 + getProb(X, output$t1, output$t0)
	}
	prob.1 <- prob.1/500

	# Combine results
	if (N ==1) {
		prob.hit <- prob.1
	}
	else {
		prob.hit <- cbind(prob.hit, prob.1)
	}
}
write.csv(prob.hit,'result/prob_hit_enzyme1.csv')


### Test random peptides' prob being hit
# Initialization
X <- data
Z <- c(1:13)
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
rpeptide <- ceiling(matrix(6*runif(500*21), ncol=21))
prob <- rep(0, 500)
for (t in 1:500) {
	output <- gibbsIter(output$Y, output$W, output$X, output$Z)
	prob <- prob + getProb(rpeptide, output$t1, output$t0)
}
prob <- prob/500
write.csv(prob, 'result/prob_random_enzyme1.csv')