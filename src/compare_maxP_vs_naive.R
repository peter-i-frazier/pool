rm(list=ls())

#=================================================================================
#Specify Paths and working directory
dataFile <- '../data/newData.csv'
classFile <- '../data/Reduced_AA_Alphabet.csv'
#=================================================================================
#import module
source('Naive_Bayes/Naive_Bayes_util.R')
source('Opt_Search/Opt_Search_util.R')

#=================================================================================
#get data 
data_org <- data.frame(read.csv(dataFile, header = T, as.is = T, sep = ","))
AAclass <- read.csv(classFile, header=T, as.is = T, sep=",")
#=================================================================================
#Set parameters
nL <- 19
nR <- 19
trainData <- getFeatures(data_org,AAclass,nL,nR)
nAA <- max(AAclass[1,])
#=================================================================================
## For AcpH
outcome <- 'AcpH'
X <- as.matrix(trainData[,c(1:(nL+nR))])
Y <- trainData[,outcome]

Gamma_0 <- 1000
Gamma_1 <- 0.05
# peptide library
nF <- dim(X)[2]
Nlib <- 1e5
maxLen <- 15
minLen <- 8
peplib <- matrix(-1, nrow=Nlib, ncol=nF)
for (n in 1:Nlib) {
	tL <- ceiling(runif(1)*(maxLen - minLen)) + minLen
	ratioL <- runif(1)*0.4 + 0.2
	L <- floor(tL * ratioL)
	R <- tL - L
	peplib[n, (nF/2-L+1):(nF/2+R)] <- ceiling(runif(L+R)*nAA)
}
maxL <- ceiling(maxLen*0.6)
maxR <- ceiling(maxLen*0.8)
#recom
recom <- maxP_search3(X, Y, AAclass, Nrec=40, itr=400, peplib, maxL=maxL, maxR=maxR, Gamma_1 = Gamma_1, Gamma_0 = Gamma_0, add_ins = 10)