library(MCMCpack)
library(Rlab)

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
# source('../Naive_Bayes/Naive_Bayes.R')

#### NOTE!!!! Feature matrix MUST NOT have NAs, for NAs, replace with negative values, like -1

sparsePrior <- function(X, Y, test, nAA, burnin.step, record.step) {
# some constant
nF <- dim(X)[2]
X.1 <- X[Y==1,]
X.0 <- X[Y==0,]
factor <- length(Y==1)/length(Y==0)
# Initialization
ztable <- Ztable(nAA)
P.1 <- rep(0.5, nF)
P.0 <- P.1
Z.1 <- c()
for (i in 1:nF) {
	Z.1 <- cbind(Z.1,rbern(nAA, P.1[i]))
}
Z.0 <- Z.1
# burn in step
for (t in 1:burnin.step) {
	mu.1 <- getMu(X.1, factor, Z.1, P.1)
	mu.0 <- getMu(X.0, 1, Z.0, P.0)
	Z.1 <- getZ(X.1, P.1, mu.1, ztable)
	Z.0 <- getZ(X.0, P.0, mu.0, ztable)
	P.1 <- getP(Z.1)
	P.0 <- getP(Z.0)
}
# record step
prob <- c()
for (t in 1:record.step) {
	mu.1 <- getMu(X.1, factor, Z.1, P.1)
	mu.0 <- getMu(X.0, 1, Z.0, P.0)
	Z.1 <- getZ(X.1, P.1, mu.1, ztable)
	Z.0 <- getZ(X.0, P.0, mu.0, ztable)
	P.1 <- getP(Z.1)
	P.0 <- getP(Z.0)
	if (is.vector(test)) {
		prob <- c(prob, getProb(test, mu.1, mu.0, Z.1, Z.0))
	} else {
		prob <- rbind(prob, getProb(test, mu.1, mu.0, Z.1, Z.0))
	}
}
return (prob)
}

getMu <- function(train, factor, Z, P) {
	nAA <- dim(Z)[1]
	nF <- dim(Z)[2]
	mu <- matrix(NA, nrow=nAA, ncol=nF)
	for (j in 1:nF) {
		distance <- as.numeric(paste(unlist(strsplit(colnames(train)[j],''))[-1],collapse=""))
		no.true <- sum(Z[,j]==1)
		if (no.true==0) {
			mu[,j] <- rep(distance**0.5, nAA) * factor
		}
		else {
			index.true <- (1:nAA)[Z[,j]==1]
			alpha <- c()
			for (i in index.true) {
				alpha <- c(alpha, distance**0.5 * factor + sum(train[,j]==i))
			}
			mu[index.true,j] <- rdirichlet(1,alpha) * rgamma(1,sum(alpha))
			mu[Z[,j]==0,j] <- rep(distance**0.5, nAA-no.true) * factor
		}
	}
	return (mu)
}

Ztable <- function(nAA){
	table <- rbind(c(0,0),c(0,1),c(1,0),c(1,1))
	for (i in 1:(nAA-2)) {
		newtable <- c()
		N <- dim(table)[1]
		for (n in 1:N) {
			newtable <- rbind(newtable, c(table[n,],0))
			newtable <- rbind(newtable, c(table[n,],1))
		}
		table <- newtable
	}
	table <- table[-1,]
	return (table)
}
getZ <- function(train, P, mu, ztable) {
	nAA <- dim(mu)[1]
	nF <- dim(mu)[2]
	Z <- matrix(NA, nrow=nAA, ncol=nF)
	N <- dim(ztable)[1]
	for (j in 1:nF) {
		prob <- rep(-1, N)
		for (n in 1:N) {
			z <- ztable[n,]
			theta <- z * mu[,j] / sum(z * mu[,j])
			# print (theta)
			# print (train[,j])
			# print (theta[train[,j]])
			prob[n] <- prod(theta[train[,j]])
		}
		# print (prob)
		prob <- prob / sum(prob)
		# print (prob)
		which.z <- 0
		pp <- runif(1)
		while (pp >0) {
			which.z <- which.z +1
			pp <- pp - prob[which.z]
		}
		Z[,j] <- ztable[which.z,]
	}
	return (Z)
}

getP <- function(Z) {
	nF <- dim(Z)[2]
	P <- rep(NA, nF)
	for (j in 1:nF) {
		P[j] <- rbeta(1, 1+ sum(Z[,j]==1), 1+ sum(Z[,j]==0))
	}
	return (P)
}

getProb <- function(test.data, mu.1, mu.0, Z.1, Z.0) {
	P <- 1e-4
	theta.1 <- c()
	theta.0 <- c()
	nF <- dim(mu.1)[2]
	for (j in 1:nF) {
		theta.1 <- cbind(theta.1, mu.1[,j]*Z.1[,j]/sum(mu.1[,j]*Z.1[,j]))
		theta.0 <- cbind(theta.0, mu.0[,j]*Z.0[,j]/sum(mu.0[,j]*Z.0[,j]))
	}
	if (is.vector(test.data)) {
			prod.1 <- 1
			prod.0 <- 1
			for (j in 1:nF) {
				if (test.data[j]!=-1) {
					prod.1 <- prod.1 * theta.1[test.data[j],j]
					prod.0 <- prod.0 * theta.0[test.data[j],j]
				}
			}
			prob <- P * prod.1/ (P * prod.1 + (1.0-P) * prod.0)
			return (prob)
		} else {
			R <- dim(test.data)[1]
			prob <- rep(0,R)
			for (r in 1:R) {
				prod.1 <- 1
				prod.0 <- 1
				for (j in 1:nF) {
					if (test.data[r,j]!=-1) {
						prod.1 <- prod.1 * theta.1[test.data[r,j],j]
						prod.0 <- prod.0 * theta.0[test.data[r,j],j]
					}
				}
				prob[r] <- P * prod.1/ (P * prod.1 + (1-P) * prod.0)
			}
		}
			return (prob)
}

getFeatures <- function(data.org, classlist, nL, nR)
{
#================================================================================
#Function: getFeatures
#
#--------------------------------------------------------------------------------
#Description:    
#    Get the feature vectors of all peptides. 
#
#--------------------------------------------------------------------------------
#Input arguments:
#    data.org
#             A matrix whose rows correspond to peptides. The matrix has n+2 
#             columns. The first n show whether each peptide works with n
#             different enzymes. The last two columns contain amino-acids at the
#             left of the serene('nterms') and right('cterms').
#    classlist
#             A matrix specifying the mapping between each amino-acids and the 
#             class it belongs to.
#    nL
#             Number of amino acids at the left of the serene that will be 
#             features.
#    nR
#             Number of amino acids at the right of the serene that will be 
#             features.
#
#--------------------------------------------------------------------------------
#Return objects:
#    feature
#             A data frame whose rows correspond to peptides and columns to 
#             features. The last columns are the outcome values(whether each 
#             peptide works with enzyme 1, 2a and 2b).
#
#--------------------------------------------------------------------------------
    #nVal: number of values each feature can take
    nVal <- length(unique(as.numeric(classlist)))
    #nOUTCOME: number of outcome values
    nOUTCOME <- dim(data.org)[2]-2
    feature <- matrix(-1, nrow=dim(data.org)[1], ncol=nL+nR+nOUTCOME)
    for (r in 1:dim(data.org)[1]) {
        sequence <- unlist(strsplit(data.org[r, 'nterm'],''))
        l.seq  <- length(sequence)
        l <- min(nL,l.seq)
        for (i in 1:l) {
		feature[r,nL+1-i] <- AAclass[1,sequence[l.seq+1-i]]
	  }
        sequence <- unlist(strsplit(data.org[r, 'cterm'],''))
        l.seq  <- length(sequence)
        c <- min(nR,l.seq)
        for (i in 1:c) {
		feature[r,nL+i] <- AAclass[1,sequence[i]]
	  }
        if( nOUTCOME != 0) {
            for (i in 1:nOUTCOME) {
                feature[r,nL+nR+i] <- data.org[r,i]
            }
        }
    }    
    #outcome.names : name of outcome values
    outcome.names <- c()
    if(nOUTCOME != 0) {
        outcome.names <- colnames(data.org)[1:nOUTCOME]
    }
    feature <- as.data.frame(feature)
    if(nOUTCOME != 0){
        colnames(feature) <- c(paste('L',nL:1,sep=""),
            paste('R',1:nR,sep=""), outcome.names)
    }
    else {
        colnames(feature) <- c(paste('L',nL:1,sep=""), 
            paste('R',1:nR,sep=""))
    }
    return( feature )
}
      
