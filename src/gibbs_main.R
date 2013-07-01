rm(list=ls())

source('gibbs_util.R')

X <- as.matrix(read.csv('FFF.csv', header=F, as.is=T))
Z <- c(1:13)
W <- X
for (i in 1:30) {
	W <- rbind(W, ceiling(6 * runif(dim(W)[2])))
}
theta.1 <- matrix(1/6,nrow=6,ncol=dim(W)[2])
theta.0 <- theta.1
Y <- sampleY(theta.1, theta.0, W, Z)

for (t in 1:1000) {
	print (t)
	theta.comb <- sampleTheta(W,Y)
	theta.1 <- theta.comb$t1
	theta.0 <- theta.comb$t0
	Y <- sampleY(theta.1, theta.0, W, Z)
	W <- sampleW(Y, theta.1, theta.0, X, Z)
	Z <- sampleZ(W, X)
}
