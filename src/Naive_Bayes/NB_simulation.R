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
source(paste(srcPath, '/getTheta_MC.R',sep=""))
source(paste(srcPath, '/NB_peptideSim.R', sep=""))
source(paste(srcPath, '/NB_peptideSim_par.R', sep=""))
source(paste(srcPath,'/NB_predict.R',sep=""))
source(paste(srcPath,'/NB_predict_par.R',sep=""))
source(paste(srcPath,'/ProbImprovement.R',sep=""))
source(paste(srcPath,'/ProbImprovement_par.R',sep=""))
source(paste(srcPath,'/ExpImprovement.R',sep=""))
source(paste(srcPath,'/ExpImprovement_par.R',sep=""))
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
nPeptide <- (300-8)/2
sim_pep <- NB_peptideSim(train = train_AcpH, classlist = AAclass, 
    outcome_name = outcome, nPep = nPeptide, isHit = 1)
sim_pep <- sim_pep[,-1]
pep_recom_sim <- paste(sim_pep[,'nterm'], 'S', sim_pep[,'cterm'], sep="")

#Mutate the existing peptides
#Here we use original data, not the table of features
data_hit <- unique(data_org[data_org$AcpH == 1, ])
nTrain <- dim(data_hit)[1]
mut_pep <- matrix("",nrow = 0, ncol = 2)
colnames(mut_pep)<-c('nterm','cterm')
while(dim(mut_pep)[1] < 300-nPeptide-8) {
	n_mut <- ceiling(4*runif(1))
	sl_pep <- data_hit[ceiling(runif(1)*nTrain),c('nterm','cterm')]
	if(runif(1)<0.5) {
	    cL <- sample(c(0:4),1)
		cR <- sample(c(0:4),1)
		if(nchar(sl_pep$nterm)+nchar(sl_pep$cterm) > 13) {
			sl_pep$nterm <- paste(unlist(strsplit(sl_pep$nterm,split=""))[-c(1:(cL+1))],collapse="")
			sl_pep$cterm <- paste(unlist(strsplit(sl_pep$cterm,split=""))[c(1:(nchar(sl_pep$cterm)-cR))],collapse="")
		}
	}
	sl_pos  <- sample.int(nchar(sl_pep$nterm)+nchar(sl_pep$cterm),n_mut)
	for(j in sl_pos) {
	    if(j <= nchar(sl_pep$nterm)){
		    sl_AA <- unlist(strsplit(sl_pep$nterm,split=""))[j] 
		} else {
			sl_AA <- unlist(strsplit(sl_pep$cterm,split=""))[j-nchar(sl_pep$nterm)]
		}
	    sl_class <- as.numeric(AAclass[sl_AA])
		new_class <- sample(c(1:8)[-sl_class],1)
		new_AA <- sample(colnames(AAclass)[AAclass == new_class],1)
		if(j <= nchar(sl_pep$nterm)){
		    tmp <- unlist(strsplit(sl_pep$nterm,split=""))
			tmp[j] <- new_AA
			sl_pep$nterm <- paste(tmp,collapse="") 
		} else {
			tmp <- unlist(strsplit(sl_pep$cterm,split=""))
			tmp[j-nchar(sl_pep$nterm)] <- new_AA
			sl_pep$cterm <- paste(tmp,collapse="") }
	}
	#check if this peptide matches previously generated ones
	is.match <- FALSE
	is.match <- !(dim(sim_pep[sim_pep[,'nterm'] == sl_pep$nterm & sim_pep[,'cterm'] == sl_pep$cterm,])[1] == 0 &&
	( dim(mut_pep[mut_pep[,'nterm'] == sl_pep$nterm & mut_pep[,'cterm'] == sl_pep$cterm,])[1] == 0 ||
	dim(mut_pep)[1] == 0))
	if(!is.match) {
	    mut_pep <- rbind(mut_pep,sl_pep)
	}
}
rownames(mut_pep) <- NULL
pep_recom_split <- as.matrix(rbind(sim_pep,mut_pep))
pep_recom_mute <- paste(mut_pep$nterm,'S',mut_pep$cterm,sep="")
peptide_recom <- c(pep_recom_sim,pep_recom_mute)

fileConn <- file("peptide_recommendation.txt")
writeLines(peptide_recom,fileConn)
close(fileConn)

#Compute the probability of simulated peptides being hit using simulation
testFeat <- getFeatures(pep_recom_split[,c('nterm','cterm')],AAclass,nL,nR)
predict_mat <- c()
for(i in 1:1000) {
    theta <- getTheta_MC(alpha = alpha, classlist = AAclass)
    predict_mat <- rbind(predict_mat, NB_predict_par(testFeat, theta))
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
exp_reduc <- ExpImprovement_par(train_AcpH, testFeat, AAclass, 100)

library(doParallel)
cl <- makeCluster(4)
registerDoParallel(cl)
system.time(NB_predict_par(testFeat, theta))[3]

#KG
org_exp_impr <- ExpImprovement(train_AcpH, kg_peptides, AAclass)

