#=================================================================================
#Specify Paths and working directory
rm(list=ls())
dataPath <- 'D:/Study/Summer2013/Peptide/data.original'
dataFile <- paste(dataPath,  '/newData.csv', sep = "")
classFile <- paste(dataPath, '/Reduced_AA_Alphabet.csv', sep = "")
srcPath <- 'D:/Study/Summer2013/Peptide/src'
wrkDir <- 'D:/Study/Summer2013/Peptide/wrkDir'
setwd(wrkDir)
#=================================================================================
#Functions might be used
source(paste(srcPath, '/Sparse_prior/sparse_prior_util_par.R', sep = ""))
source(paste(srcPath, '/NaiveBayes/Naive_Bayes_util.R',sep=""))
source(paste(srcPath, '/Opt_Search/Opt_Search_util.R',sep=""))
#=================================================================================
#get data 
data_org <- data.frame(read.csv(dataFile, header = T, as.is = T, sep = ","))
AAclass <- read.csv(classFile, header=T, as.is = T, sep=",")
#=================================================================================
#Parallel Computation Packages
require(multicore)
require(snow)
require(doParallel)
require(MCMCpack)
require(Rlab)
#=================================================================================
#Parallel Computation Setting
cl <- makeCluster(4)
registerDoParallel(cl)
#=================================================================================
#Set parameters
nL <- 19
nR <- 19
trainData <- getFeatures(data_org,AAclass,nL,nR)
nAA <- max(AAclass[1,])
burnin.step <- 5000
record.step <- 10000
#=================================================================================
## For AcpH
outcome <- 'AcpH'
X <- as.matrix(trainData[,c(1:(nL+nR))])
Y <- trainData[,outcome]
tX <- X[Y==0,][-c(1:6),]

Sample <- sample.int(dim(tX)[1], size = dim(X[Y==1,])[1]-6)
trainX <- rbind(X[c(1:16),], X[Y==1,][-c(1:10),], tX[Sample,])
trainY <- c(rep(1,10), rep(0,6), rep(1,dim(X[Y==1,])[1]-10),rep(0,dim(X[Y==1,])[1]-6))



## cross validation
registerDoMC()
system.time(
prob <- foreach (n = 1:dim(X)[1], .init = c(), .combine = rbind, .packages = c('MCMCpack','Rlab'), .inorder = TRUE, .errorhandling = 'pass') %dopar% {
	sparsePrior(X[-n,], Y[-n], X[n,], nAA, burnin.step, record.step)
})[3]

prob <- c()
system.time(
for ( n in 1:dim(X)[1] ) {
	prob <- rbind(prob, sparsePrior(X[-n,], Y[-n], X[n,], nAA, burnin.step, record.step))
})[3]

rowNames(prob) <- c()
prob_pep <- rowMeans(prob)

#Plot the ROC curve
FPR <- rep(-1, 101)
TPR <- rep(-1, 101)
for( i in 1:101 ) {
	threshold <- (i-1)/100
	label <- rep(0, length(prob_pep))
	for( j in 1:length(prob_pep) ){
		if(prob_pep[j] > threshold) {
			label[j] <- 1 }
	}
	FPR[i] <- sum((Y==0)&(label==1))/sum(Y==0)
	TPR[i] <- sum((Y==1)&(label==1))/sum(Y==1)
}
TP <- sum((Y==1)&(label==1))
FN <- sum((Y==1)&(label==0))
TN <- sum((Y==0)&(label==0))
FP <- sum((Y==0)&(label==1))



jpeg('ROC_SP.jpg')
plot(x = FPR, y = TPR, type = 'l', xlim = c(0,1), ylim = c(0,1), xlab = "false positive rate", ylab = "true positive rate")
dev.off()