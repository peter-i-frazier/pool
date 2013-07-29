rm(list=ls())
source('gibbs_util.R')
source('../Naive_Bayes/getFeatures.R')
#=================================================================================
#Specify Paths and working directory
dataFile <- '../../data/binaryData_v2.csv'
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
itr <- 1000
Nrec <- 1e6
#=================================================================================
## For AcpH
outcome <- 'AcpH'
X.train <- as.matrix(trainData[,c(1:(nL+nR))])
Y.train <- trainData[,outcome]
Npos <- sum(Y.train)
Nneg <- length(Y.train) - Npos

#calculate mean and std of theta
storage.1 <- matrix(NA, nrow=(nL+nR)*nVal, ncol=itr+2)
storage.0 <- matrix(NA, nrow=(nL+nR)*nVal, ncol=itr+2)
for (n in 1:itr) {
	theta.comb <- sampleTheta(X.train, Y.train, nVal, Npos, Nneg)
	theta.1 <- theta.comb$theta_1
	theta.0 <- theta.comb$theta_0
	for (i in 1:nVal) {
		for (j in 1:(nL+nR)) {
			storage.1[(i-1)*(nL+nR)+j, n] <- theta.1[i,j]
			storage.0[(i-1)*(nL+nR)+j, n] <- theta.0[i,j]
		}
	}
}
for (i in 1:dim(storage.1)[1]) {
	storage.1[i,itr+1] <- mean(storage.1[i, 1:itr])
	storage.1[i,itr+2] <- sd(storage.1[i, 1:itr])
	storage.0[i,itr+1] <- mean(storage.0[i, 1:itr])
	storage.0[i,itr+2] <- sd(storage.0[i, 1:itr])
}
theta.1.sd <- matrix(NA, nrow=nVal, ncol=(nL+nR))
theta.0.sd <- matrix(NA, nrow=nVal, ncol=(nL+nR))
for (i in 1:dim(storage.1)[1]) {
	theta.1[ceiling(i/(nL+nR)), i-(nL+nR)*floor(i/(nL+nR))] <- storage.1[i, itr+1]
	theta.1.sd[ceiling(i/(nL+nR)), i-(nL+nR)*floor(i/(nL+nR))] <- storage.1[i, itr+2]
	theta.0[ceiling(i/(nL+nR)), i-(nL+nR)*floor(i/(nL+nR))] <- storage.0[i, itr+1]
	theta.0.sd[ceiling(i/(nL+nR)), i-(nL+nR)*floor(i/(nL+nR))] <- storage.0[i, itr+2]
}

print ('simulate recommendations')
#simulate recommendations
random.peptides <- matrix(NA, nrow=Nrec, ncol=(nL+nR))
random.nL <- ceiling(runif(Nrec) * 3) + 1    # 2-4 uniform
random.nR <- ceiling(runif(Nrec) * 7) + 2    # 3-9 uniform
for (n in 1:Nrec) {
	pp <- runif(random.nL[n]+random.nR[n])
	theta <- theta.1[,(nL-random.nL[n]+1):(nL+random.nR[n])]
	one.peptide <- rep(1,random.nL[n]+random.nR[n])
	for (j in 1:dim(theta)[1]) {
		one.peptide <- one.peptide + ((pp-theta[j,])>0)
		pp <- pp-theta[j,] 
	}
	random.peptides[n,(nL-random.nL[n]+1):(nL+random.nR[n])] <- one.peptide
}
prob.random.peptides <- getProb(random.peptides, theta.1, theta.0)
table.forSort <- cbind(prob.random.peptides, random.peptides)
sorted.table <- table.forSort[order(table.forSort[,1],decreasing=T),]
recommend.list <- sorted.table[1:121,-1]
recommend.prob <- sorted.table[1:121,1]

