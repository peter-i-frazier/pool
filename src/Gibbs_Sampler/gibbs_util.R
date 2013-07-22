library(MCMCpack)
#================================================================================
#Module: Gibbs Sampler supporting functions
#
#--------------------------------------------------------------------------------
#Description:
#    Sample from conditional distributions of theta, Y, W, Z iteratively, to get 
#	 the ture distribution of theta given data. Full algorithm see write up.
#
#================================================================================


#--------------------------------------------------------------------------------
# Calculate theta given training data
# Input: training data   Type: matrix
# Output: theta          Type: matrix
#--------------------------------------------------------------------------------
getTheta <- function(train) {
	K <- dim(train)[2]
	R <- dim(train)[1]
	theta <- c()
	for (col in 1:K) {
		count <- rep(0,6)
		for (r in 1:R) {
			count[train[r,col]] <- count[train[r,col]] + 1
			alpha <- count + rep(abs(col-12.5)**0.5,6)
		}
			theta <- cbind(theta, c(rdirichlet(1,alpha)))
	}
	return (theta)
}
sampleTheta <- function(W,Y) {
	train.1 <- W[Y==1,]
	train.0 <- W[Y==0,]
	theta.1 <- getTheta(train.1)
	theta.0 <- getTheta(train.0)
	theta.comb <- list(t1=theta.1, t0=theta.0)
	return (theta.comb)
}

sampleY <- function(theta.1, theta.0, W, Z) {
	P <- 1e-4
	M <- dim(W)[1]
	K <- dim(W)[2]
	Y <- rep(0,M)
	for (m in 1:M) {
		if (sum(m==Z)>0) {
			Y[m] <- 1
		}
		else {
			x <- W[m,]
			prod.1 <- 1
			prod.0 <- 1
			for (j in 1:K) {
				prod.1 <- prod.1 * theta.1[x[j],j]
				prod.0 <- prod.0 * theta.0[x[j],j]
			}
			prob <- P * prod.1/ (P * prod.1 + (1.-P) * prod.0)
			if (runif(1) < prob) {Y[m] <- 1}
			else {Y[m] <- 0}
		}
	}
	return (Y)
}

sampleW <- function(Y, theta.1, theta.0, data, Z) {
	M <- length(Y)
	K <- dim(data)[2]
	N <- dim(data)[1]
	W <- matrix(1, nrow=M, ncol=K)
	for (m in 1:M) {
		sample.or.not <- 1
		for (n in 1:N) {
			if (Z[n]==m) {
				W[m,] <- data[n,]
				sample.or.not <- 0
				break
			}
		}
		if (sample.or.not) {
			W[m,] <- rep(1,K)
			if (Y[m]==1) {
				pp <- runif(K)
				for (j in 1:dim(theta.1)[1]) {
					W[m,] = W[m,] + ((pp-theta.1[j,])>0)
					pp <- pp-theta.1[j,]
				}
			}
			else {
				pp <- runif(K)
				for (j in 1:dim(theta.0)[1]) {
					W[m,] = W[m,] + ((pp-theta.0[j,])>0)
					pp <- pp-theta.0[j,]
				}
			}
		}
	}
	return (W)
}

sampleZ <- function(W, data) {
	M <- dim(W)[1]
	N <- dim(data)[1]
	Z <- rep(0,N)
	for (n in 1:N) {
		candidates <- c()
		for (m in 1:M) {
			if (prod(data[n,]==W[m,])) {candidates <- c(candidates, m)}
		}
		Len <- length(candidates)
		Z[n] <- candidates[ceiling(Len * runif(1))]
	}
	return (Z)
}

getProb <- function(X, theta.1, theta.0) {
		P <- 1e-4
		K <- dim(X)[2]
		prob <- rep(0,dim(X)[1])
		for (i in 1:dim(X)[1]) {
			x <- X[i,]
			prod.1 <- 1
			prod.0 <- 1
			for (j in 1:K) {
				prod.1 <- prod.1 * theta.1[x[j],j]
				prod.0 <- prod.0 * theta.0[x[j],j]
			}
			prob[i] <- P * prod.1/ (P * prod.1 + (1.-P) * prod.0)
		}
		return (prob)
}

gibbsIter <- function(Y, W, X, Z, itr=1) {
	if (itr == 1) {
		theta.comb <- sampleTheta(W,Y)
		theta.1 <- theta.comb$t1
		theta.0 <- theta.comb$t0
		Y <- sampleY(theta.1, theta.0, W, Z)
		W <- sampleW(Y, theta.1, theta.0, X, Z)
		Z <- sampleZ(W, X)
		result <- list(Y=Y, W=W, X=X, Z=Z, t1=theta.1, t0=theta.0)
		return (result)
	}
	else {
		for (n in 1:itr) {
			theta.comb <- sampleTheta(W,Y)
			theta.1 <- theta.comb$t1
			theta.0 <- theta.comb$t0
			Y <- sampleY(theta.1, theta.0, W, Z)
			W <- sampleW(Y, theta.1, theta.0, X, Z)
			Z <- sampleZ(W, X)
		}
		result <- list(Y=Y, W=W, X=X, Z=Z, t1=theta.1, t0=theta.0)
		return (result)
	}
}