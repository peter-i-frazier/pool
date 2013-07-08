#=================================================================================
#Specify Paths and working directory
rm(list=ls())
dataPath <- 'D:/Study/Summer2013/Peptide/data.original'
dataFile <- paste(dataPath,  '/binaryData_v2.csv', sep = "")
classFile <- paste(dataPath, '/Reduced_AA_Alphabet.csv', sep = "")
srcPath <- 'D:/Study/Summer2013/Peptide/src'
wrkDir <- 'D:/Study/Summer2013/Peptide/wrkDir'
setwd(wrkDir)
#=================================================================================
#Functions might be used
source(paste(srcPath, '/getFeatures.R',sep=""))
source(paste(srcPath, '/Dirichlet_Parameter.R',sep=""))
source(paste(srcPath, '/getTheta_MC.R',sep=""))
source(paste(srcPath, '/rdist.R', sep=""))
source(paste(srcPath, '/NB_peptideSim.R', sep=""))
source(paste(srcPath,'/NB_predict.R',sep=""))
source(paste(srcPath,'/ProbImprovement.R',sep=""))
source(paste(srcPath,'/ExpImprovement.R',sep=""))
#=================================================================================
#get data 
data <- data.frame(read.csv(dataFile, header = T, as.is = T, sep = ","))
class <- read.csv(classFile, header=T, as.is = T, sep=",")
#=================================================================================
#Set parameters
nL <- 19
nR <- 19

trainData <- getFeatures(data,class,nL,nR)
nVal <- max(class)
outcome <- 'AcpH'
train_AcpH <- trainData[,c(1:(nL+nR),which(colnames(trainData)==outcome))]
# Get the alpha parameter of posterior Dirichlet distribution
alpha <- Dirichlet_Parameter(train_AcpH, nVal)
#=================================================================================
#Simulate peptides and write them to a text file
nPeptide <- 300
pep_recom_split <- NB_peptideSim(train = train_AcpH, nVal = nVal, outcome = outcome, nPep = nPeptide, isHit = 1)
peptide_recom <- paste(pep_recom_split[,2], 'S', pep_recom_split[,3], sep="")

fileConn <- file("peptide_recommendation.txt")
writeLines(peptide_recom,fileConn)
close(fileConn)

#Compute the probability of simulated peptides being hit using simulation
testFeat <- getFeatures(pep_recom_split[,c('nterm','cterm')],class,nL,nR)
predict_mat <- c()
for(i in 1:1000) {
    theta <- getTheta_MC(alpha = alpha, nVal = nVal)
    predict_mat <- rbind(predict_mat,
        NB_predict(testFeat, theta, prior.positive = 0.5))
}
predict <- colMeans(predict_mat)

#Compute probability of improvement using our simulate peptide_recom
library(Rlab)
prob.improvement <- ProbImprovement(train_AcpH, testFeat, 100, nVal, prior.positive = 0.5)





