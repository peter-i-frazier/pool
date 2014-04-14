#=================================================================================
#Specify Paths and working directory
rm(list=ls())
dataPath <- '/fs/home/py75/Documents/peptide/data'
dataFile <- paste(dataPath,  '/binaryData_v2.csv', sep = "")
classFile <- paste(dataPath, '/Reduced_AA_Alphabet.csv', sep = "")
srcPath <- '/fs/home/py75/Documents/peptide/src'
wrkDir <- '/fs/home/py75/Documents/peptide/wrkDir'
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
#Parallel Computation Packages
require(multicore)
require(snow)
require(doParallel)
#=================================================================================
#Parallel Computation Setting
hosts <- c(rep("whale",52),rep("ahab",24),rep("flask",24),rep("starbuck",24))
cl <- makeCluster(hosts, type = "SOCK")
registerDoParallel(cl)
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
nPeptide <- (250-8)/2
sim_pep <- NB_peptideSim(train = train_AcpH, classlist = AAclass, 
    outcome_name = outcome, nPep = nPeptide, isHit = 1)
sim_pep <- sim_pep[,c('nterm','cterm')]
#Mutate the existing peptides
#Here we use original data, not the table of features
data_hit <- unique(data_org[data_org$AcpH == 1, ])
nTrain <- dim(data_hit)[1]
mut_pep <- matrix("",nrow = 0, ncol = 2)
colnames(mut_pep)<-c('nterm','cterm')
while(dim(mut_pep)[1] < 250-nPeptide-8) {
	n_mut <- sample.int(4,1)
	
	sl_pep <- data_hit[ceiling(runif(1)*nTrain),c('nterm','cterm')]
	cL <- sample(c(6:9),1)
	cR <- 19 - cL
	if(nchar(sl_pep$nterm)+nchar(sl_pep$cterm) > 13) {
		sl_pep$nterm <- paste(unlist(strsplit(sl_pep$nterm,split=""))[c((nchar(sl_pep$nterm)-cL+1):nchar(sl_pep$nterm))],collapse="")
		sl_pep$cterm <- paste(unlist(strsplit(sl_pep$cterm,split=""))[c(1:cR)],collapse="")
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
#==================================================================================
#Compute the probability of simulated peptides being hit using simulation
testFeat <- getFeatures(pep_recom_split[,c('nterm','cterm')],AAclass,nL,nR)
predict_mat <- c()
predict_mat <- foreach(i=1:1000, .combine = rbind, .multicombine = TRUE, .init = c(),
       .export = c("getTheta_MC","NB_predict")) %dopar% {
    theta <- getTheta_MC(alpha = alpha, classlist = AAclass)
    NB_predict(testFeat, theta)
}
predict <- colMeans(predict_mat)
#Compute probability of improvement using our simulate peptide_recom
prob_impr <- ProbImprovement_par(train_AcpH, testFeat, AAclass, 1000)
#Compute expected reduction in length using simulated peptide_recom
exp_reduc <- ExpImprovement_par(train_AcpH, testFeat, AAclass, 1000)

