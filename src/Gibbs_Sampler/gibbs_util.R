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
source('../Naive_Bayes/Naive_Bayes.R')

getTheta <- function(train, nVal) {
	K <- dim(train)[2]
	R <- dim(train)[1]
	theta <- c()
	for (col in 1:K) {
		count <- rep(0,nVal)
		distance <- as.numeric(paste(unlist(strsplit(colnames(train)[col],''))[-1],collapse=""))
		for (r in 1:R) {
			if (!is.na(train[r,col])) {
				count[train[r,col]] <- count[train[r,col]] + 1
			}
		}
		alpha <- count + rep(distance**0.5,nVal)
		theta <- cbind(theta, c(rdirichlet(1,alpha)))
	}
	return (theta)
}
sampleTheta <- function(W, Y, nVal) {
	train.1 <- W[Y==1,]
	train.0 <- W[Y==0,]
	theta.1 <- getTheta(train.1, nVal)
	theta.0 <- getTheta(train.0, nVal)
	theta.comb <- list(theta_1=theta.1, theta_0=theta.0)
	return (theta.comb)
}

# sampleTheta <- function(W,Y, AAclass) {
# 	alpha <- Dirichlet_Parameter(cbind(W,Y),AAclass)
# 	theta.comb <- getTheta_MC(cbind(W,Y), alpha, AAclass)
# 	return (theta.comb)
# }

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
			prod.1 <- 1
			prod.0 <- 1
			for (j in 1:K) {
				if (!is.na(W[m,j])) {
					prod.1 <- prod.1 * theta.1[W[m,j],j]
					prod.0 <- prod.0 * theta.0[W[m,j],j]
				}
			}
			prob <- P * prod.1/ (P * prod.1 + (1-P) * prod.0)
			if (runif(1) < prob) {Y[m] <- 1}
			else {Y[m] <- 0}
		}
	}
	return (Y)
}

sampleW <- function(Y, theta.1, theta.0, Data, Z) {
	M <- length(Y)
	K <- dim(Data)[2]
	N <- dim(Data)[1]
	W <- matrix(NA, nrow=M, ncol=K)
	colnames(W) <- colnames(Data)
	for (m in 1:M) {
		sample.or.not <- 1
		for (n in 1:N) {
			if (Z[n]==m) {
				W[m,] <- Data[n,]
				sample.or.not <- 0
				break
			}
		}
		if (sample.or.not) {
			W[m,] <- rep(1,K)
			if (Y[m]==1) {
				pp <- runif(K)
				for (j in 1:dim(theta.1)[1]) {
					W[m,] <- W[m,] + ((pp-theta.1[j,])>0)
					pp <- pp-theta.1[j,]
				}
			}
			else {
				pp <- runif(K)
				for (j in 1:dim(theta.0)[1]) {
					W[m,] <- W[m,] + ((pp-theta.0[j,])>0)
					pp <- pp-theta.0[j,]
				}
			}
		}
	}
	return (W)
}

sampleZ <- function(W, Data) {
	M <- dim(W)[1]
	N <- dim(Data)[1]
	Z <- rep(0,N)
	for (n in 1:N) {
		candidates <- c()
		for (m in 1:M) {
			isnot.same <- 0
			diff <- W[m,] - Data[n,]
			for (j in 1:dim(W)[2]) {
				if (!is.na(diff[j])) {
					isnot.same <- isnot.same + abs(diff[j])
				}
			}
			if (isnot.same==0) {
				candidates <- c(candidates, m)
			}
		}
		Len <- length(candidates)
		Z[n] <- candidates[ceiling(Len * runif(1))]
	}
	return (Z)
}

getProb <- function(test.data, theta.1, theta.0) {
		P <- 1e-4
		if (is.vector(test.data)) {
			K <- length(test.data)
			prod.1 <- 1
			prod.0 <- 1
			for (j in 1:K) {
				if (!is.na(test.data[j])) {
					prod.1 <- prod.1 * theta.1[test.data[j],j]
					prod.0 <- prod.0 * theta.0[test.data[j],j]
				}
			}
			prob <- P * prod.1/ (P * prod.1 + (1.-P) * prod.0)
			return (prob)
		} else {
			K <- dim(test.data)[2]
			R <- dim(test.data)[1]
			prob <- rep(0,R)
			for (r in 1:R) {
				prod.1 <- 1
				prod.0 <- 1
				for (j in 1:K) {
					if (!is.na(test.data[r,j])) {
						prod.1 <- prod.1 * theta.1[test.data[r,j],j]
						prod.0 <- prod.0 * theta.0[test.data[r,j],j]
					}
				}
				prob[r] <- P * prod.1/ (P * prod.1 + (1-P) * prod.0)
			}
			return (prob)
		}
}

gibbsSampler <- function(train.data, test.data, nVal, M, burnin.step, record.step) {
	# some constant
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
		theta.comb <- sampleTheta(W, Y, nVal)
		theta.1 <- theta.comb$theta_1
		theta.0 <- theta.comb$theta_0
		Y <- sampleY(theta.1, theta.0, W, Z)
		W <- sampleW(Y, theta.1, theta.0, train.data, Z)
		Z <- sampleZ(W, train.data)
		prob <- prob + getProb(test.data, theta.1, theta.0)
	}
	prob <- prob/record.step
	return (prob)
}

