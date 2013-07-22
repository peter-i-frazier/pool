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
source(paste(srcPath, '/NB_peptideSim_par.R', sep=""))
source(paste(srcPath,'/NB_predict.R',sep=""))
source(paste(srcPath,'/NB_predict_par.R',sep=""))
source(paste(srcPath,'/ProbImprovement.R',sep=""))
source(paste(srcPath,'/ProbImprovement_par.R',sep=""))
source(paste(srcPath,'/ExpImprovement.R',sep=""))
#=================================================================================
#get data 
data_org <- data.frame(read.csv(dataFile, header = T, as.is = T, sep = ","))
AAclass <- read.csv(classFile, header=T, as.is = T, sep=",")
#=================================================================================
#Set parameters
nL <- 19
nR <- 19

trainData <- getFeatures(data_org,AAclass,nL,nR)
outcome <- 'AcpH'
train_AcpH <- trainData[,c(1:(nL+nR),which(colnames(trainData)==outcome))]
# Get the alpha parameter of posterior Dirichlet distribution
alpha <- Dirichlet_Parameter(train_AcpH, AAclass)
#=================================================================================
#Simulate peptides and write them to a text file
nPeptide <- 300
pep_recom_split <- NB_peptideSim(train = train_AcpH, classlist = AAclass, 
    outcome_name = outcome, nPep = nPeptide, isHit = 1)
peptide_recom <- paste(pep_recom_split[,2], 'S', pep_recom_split[,3], sep="")

fileConn <- file("peptide_recommendation.txt")
writeLines(peptide_recom,fileConn)
close(fileConn)

#Compute the probability of simulated peptides being hit using simulation
testFeat <- getFeatures(pep_recom_split[,c('nterm','cterm')],AAclass,nL,nR)
predict_mat <- c()
for(i in 1:1000) {
    theta <- getTheta_MC(alpha = alpha, classlist = AAclass)
    predict_mat <- rbind(predict_mat,
        NB_predict(testFeat, theta,prior.positive = 10**-3))
}
predict <- colMeans(predict_mat)

#Compute probability of improvement using our simulate peptide_recom
#Bernoulli random number generating function rbern(...) requires the 
#package Rlab
library(Rlab)
prob_impr <- ProbImprovement(train_AcpH, testFeat, AAclass, 100)
prob_impr_par <- ProbImprovement_par(train_AcpH, testFeat, AAclass, 100)
#Compute expected reduction in length using simulated peptide_recom
exp_reduc <- ExpImprovement(train_AcpH, testFeat, AAclass, 100)

library(doParallel)
cl <- makeCluster(4)
registerDoParallel(cl)
system.time(test_Pep <- NB_peptideSim_par(train = train_AcpH, outcome_name = 'AcpH', classlist = AAclass, nPep = 10000))[3]

#KG
org_exp_impr <- ExpImprovement(train_AcpH, kg_peptides, AAclass)

