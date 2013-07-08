library(gdata)
library(MCMCpack)
rm(list=ls())

# AAclass <- read.xls('../data/AAclass.xlsx', header=T, as.is=T)
# binary.data <- data.frame(read.xls('../data/binaryData_up.xlsx', header=T, as.is=T))
# ez2.feature <- matrix(0, nrow=dim(binary.data[1]), ncol=38)
# for (r in 1:dim(binary.data)[1]) {
# 	sequence <- unlist(strsplit(binary.data[r, 'nterm'],''))
# 	for (i in 1:length(sequence)) {
# 		ez2.feature[r,i] <- AAclass[1,sequence[i]]
# 	}
# 	sequence <- unlist(strsplit(binary.data[r, 'cterm'],''))
# 	for (i in 1:length(sequence)) {
# 		ez2.feature[r,i+19] <- AAclass[1,sequence[i]]
# 	}
# }
# ez2.1.x <- ez2.feature[1:8,]
# ez2.0.x <- ez2.feature[9:16,]

Feature <- read.csv('feature.csv', header=F, as.is=T)
ez2.feature <- as.matrix(Feature[-1,])
ez2.1.x <- ez2.feature[1:8,]
ez2.0.x <- ez2.feature[9:16,]

source('getTheta.R')
#theta with full data
theta.0 <- getTheta(ez2.0.x)
theta.1 <- getTheta(ez2.1.x)
write.csv(theta.1/theta.0, 'ratio_theta_reduced.csv')

#cross validation
getAlpha.matrix <- function(train) {
	K <- dim(train)[2]
	R <- dim(train)[1]
	alpha.matrix <- c()

	for (col in 1:K) {
		count <- rep(0,6)
		for (r in 1:R) {
			count[train[r,col]] <- count[train[r,col]] + 1
			alpha <- count + rep(abs(Feature[1,col]-19.5)**0.25,6)
		}
			alpha.matrix <- cbind(alpha.matrix, alpha)
	}
	return (alpha.matrix)
}
#alpha with full data
alpha.0 <- getAlpha.matrix(ez2.0.x)
alpha.1 <- getAlpha.matrix(ez2.1.x)
N <- 10000
mytest.prob <- matrix(0, nrow=dim(ez2.feature)[1], ncol=dim(ez2.feature)[2])
#for y=1
test.prob.1 <- rep(1,dim(ez2.1.x)[1])
for (n in 1:dim(ez2.1.x)[1]) {
	train <- ez2.1.x[-n,]
	compensate.row <- ceiling(dim(train)[1] * runif(1))
	train <- rbind(train[compensate.row,], train)
	# train <- ez2.1.x
	test <- ez2.1.x[n,]
	alpha.matrix.1 <- getAlpha.matrix(train)
	alpha.matrix.0 <- alpha.0
	theta.sample.1 <- c()
	theta.sample.0 <- c()
	for (i in 1:length(test)) {
		theta.sample.1 <- cbind(theta.sample.1, rdirichlet(N,alpha.matrix.1[,i])[,test[i]])
		theta.sample.0 <- cbind(theta.sample.0, rdirichlet(N,alpha.matrix.0[,i])[,test[i]])
	}
	prob.theta <- rep(0,N)
	for (i in 1:N) {
		prob.theta[i] <- prod(theta.sample.1[i,])/(prod(theta.sample.1[i,])+prod(theta.sample.0[i,]))
	}
	test.prob.1[n] <- mean(prob.theta)
	for (i in 1:length(test)) {
		mytest.prob[n,i] <- mean(theta.sample.1[,i]/theta.sample.0[,i])
	}
}

#for y=0
test.prob.0 <- rep(1,dim(ez2.0.x)[1])
for (n in 1:dim(ez2.0.x)[1]) {
	train <- ez2.0.x[-n,]
	compensate.row <- ceiling(dim(train)[1] * runif(1))
	train <- rbind(train[compensate.row,], train)
	# train <- ez2.0.x
	test <- ez2.0.x[n,]
	alpha.matrix.0 <- getAlpha.matrix(train)
	alpha.matrix.1 <- alpha.1
	theta.sample.1 <- c()
	theta.sample.0 <- c()
	for (i in 1:length(test)) {
		theta.sample.1 <- cbind(theta.sample.1, rdirichlet(N,alpha.matrix.1[,i])[,test[i]])
		theta.sample.0 <- cbind(theta.sample.0, rdirichlet(N,alpha.matrix.0[,i])[,test[i]])
	}
	prob.theta <- rep(0,N)
	for (i in 1:N) {
		prob.theta[i] <- prod(theta.sample.1[i,])/(prod(theta.sample.1[i,])+prod(theta.sample.0[i,]))
	}
	test.prob.0[n] <- mean(prob.theta)
	for (i in 1:length(test)) {
		mytest.prob[n+dim(ez2.1.x)[1],i] <- mean(theta.sample.1[,i]/theta.sample.0[,i])
	}
}
write.csv(mytest.prob,'mytest_prob_of_eachAA.csv')
test.prob <- c(test.prob.1,test.prob.0)
write.csv(test.prob, 'prob_CV_reducedAA.csv')