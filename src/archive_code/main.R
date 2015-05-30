#######			Test here

#=================================================================================
#Specify Paths and working directory
rm(list=ls())
dataPath <- 'D:/Study/Summer2013/Peptide/data'
dataFile <- paste(dataPath,  '/Data_07.csv', sep = "")
classFile <- paste(dataPath, '/Reduced_AA_Alphabet.csv', sep = "")
srcPath <- 'D:/Study/Summer2013/Peptide/src'
wrkDir <- 'D:/Study/Summer2013/Peptide/wrkDir'
setwd(wrkDir)
#=================================================================================
#Functions might be used
source(paste(srcPath, '/SparsePrior/sparse_prior_util_par.R', sep = ""))
source(paste(srcPath, '/NaiveBayes/Naive_Bayes_util.R',sep=""))
source(paste(srcPath, '/OptSearch/Opt_Search_util.R',sep=""))
source(paste(srcPath, '/writePep.R', sep = ""))


#=================================================================================


#######			Run on Cluster

#=================================================================================
#Specify Paths and working directory
rm(list=ls())
dataPath <- '/fs/home/py75/Documents/peptide/data'
dataFile <- paste(dataPath,  '/Data_10.csv', sep = "")
classFile <- paste(dataPath, '/Reduced_AA_Alphabet.csv', sep = "")
srcPath <- '/fs/home/py75/Documents/peptide/src'
wrkDir <- '/fs/home/py75/Documents/peptide/wrkDir'
setwd(wrkDir)
#=================================================================================
#Functions might be used
source(paste(srcPath, '/SparsePrior/old_sparse_prior_util.R', sep = ""))
source(paste(srcPath, '/NaiveBayes/Naive_Bayes_util.R',sep=""))
source(paste(srcPath, '/OptSearch/Opt_Search_util.R',sep=""))
source(paste(srcPath, '/writePep.R', sep = ""))
source(paste(srcPath, '/ROC_plot.R', sep = ""))
#=================================================================================
#Parallel Computation Packages
require(multicore)
require(snow)
require(doParallel)
require(MCMCpack)
require(Rlab)
#=================================================================================
#Parallel Computation Setting
hosts <- rep('whale',40)
hosts <- c(rep("whale",52),rep("ahab",24),rep("flask",24),rep("starbuck",24))
hosts <- c(rep("ishmael",28),rep("stubb",30),rep("daggoo",30),rep("tashtego",30))
cl <- makeCluster(hosts, type = "SOCK")

cl <- makeCluster(4)
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
nAA <- max(AAclass[1,])
#=================================================================================
## For AcpH
outcome <- 'AcpH'
X <- as.matrix(trainData[,c(1:(nL+nR))])
Y <- trainData[,outcome]


#########		Cross Validation
##Vanilla Naive Bayes
prob_vanilla <- foreach (n = 1:dim(X)[1], .init = c(), .combine = c, .packages = c('MCMCpack','Rlab'), .inorder = TRUE, .errorhandling = 'pass') %dopar% {
	Naive_Bayes(X[-n,], Y[-n], X[n,], AAclass, Gamma_0 = 1000, Gamma_1 = 0.05, predIter = 500)
}


##Sparse Prior
prob_sparse <- foreach (n = 1:dim(X)[1], .init = c(), .combine = rbind, .packages = c('MCMCpack','Rlab'), .inorder = TRUE, .errorhandling = 'pass') %dopar% {
	sparsePrior(X[-n,], Y[-n], X[n,], nAA, burnin.step, record.step)
}


#########		Find optimal Gamma parameter pair (Gamma_0, Gamma_1) using cross validation		##########
Gamma_1 <- c(0.0001, 0.001, 0.01, 0.05, 0.1, 0.25, 0.5)
Gamma_0 <- c(10**c(1:6))
auc <- matrix(-1, nrow = length(Gamma_1), ncol = length(Gamma_0))
for(i in 1:length(Gamma_1)){
	for(j in 1:length(Gamma_0)){
		prob_vanilla <- foreach (n = 1:dim(X)[1], .init = c(), .combine = c, .packages = c('MCMCpack','Rlab'), .inorder = TRUE, .errorhandling = 'pass') %dopar% {
			Naive_Bayes(X[-n,], Y[-n], X[n,], AAclass, Gamma_0 = Gamma_0[j], Gamma_1 =Gamma_1[i])
		}
		FPR <- rep(-1, length(prob_vanilla))
		TPR <- rep(-1, length(prob_vanilla))
		thresholds <- sort(prob_vanilla)
		for( k in 1:length(prob_vanilla) ) {
				threshold <- thresholds[k]
				label <- rep(0, length(prob_vanilla))
				for( l in 1:length(prob_vanilla) ){
					if(prob_vanilla[l] >= threshold) {
						label[l] <- 1 }
				}
		FPR[k] <- sum((Y==0)&(label==1))/sum(Y==0)
		TPR[k] <- sum((Y==1)&(label==1))/sum(Y==1)
		}	
	auc[i,j] = AUC(rev(FPR), rev(TPR))
	}
}

Gamma_0 = 1000
Gamma_1 = 0.05
#recom_1: Simulation. Add 10 each time. Nlib = 1e5
recom_3 <- maxP_search(X, Y, AAclass, Nrec=10, itr=1000, Nlib = 1e5, maxLen = 10, minLen = 5, Gamma_0 = Gamma_0, Gamma_1 = Gamma_1, add_ins = 10)

#recom_1_prime: Simulation. Add 10 each time. Nlib = 1e6
recom_1_prime <- maxP_search(X, Y, AAclass, Nrec=40, itr=400, Nlib = 1e6, maxL=9, maxR=10, add_ins = 10, rootpath=NA)

#recom_2: Simulation. Add one each time
recom_4 <- maxP_search(X, Y, AAclass, Nrec=40, itr=400, Nlib = 1e5, maxLen = 10, minLen = 5, Gamma_0 = Gamma_0, Gamma_1 = Gamma_1, add_ins = 1)

#recom_3: Mean. Add 10 each time
recom_3 <- maxP_search_2(X, Y, AAclass, Nrec = 10, maxLen = 10, minLen = 5, Gamma_0 = Gamma_0, Gamma_1 = Gamma_1, add_ins = 20)

#recom_4: Mean. Add 1 each time
recom_4 <- maxP_search_2(X, Y, AAclass, Nrec = 58, maxLen = 10, minLen = 5, Gamma_0 = Gamma_0, Gamma_1 = Gamma_1, add_ins = 1)

#recom_5: Mutation

recom_all <- rbind(recom_1, recom_2, recom_3, recom_4)
writePep(X, Y, recom_all, AAclass, Gamma_0, Gamma_1)

#Plot Expected Improvement vs. n for each recommendation set.
exImp_1 <- c()
for (i in 2:dim(recom_1)[1]){
	exImp_1 <- c(exImp_1, ExpImprovement_par(X, Y, recom_1[1:i,], AAclass)$improve)
}
jpeg('exImp_1.jpg')
plot(x = c(2:dim(recom_1)[1]), y = exImp_1, type = 'l', xlab = "number peptides included", ylab = "expected improvement", main = "Simulated theta, 10 repeats each new recom")
dev.off()


#predict probability of being hit
alpha <- Dirichlet_Parameter(X, Y, AAclass, Gamma_0 = 1.5, Gamma_1 = 0.05)
pred_1 <- foreach(i=1:1000, .init = c(), .combine = rbind) %dopar% {
	theta <- getTheta_MC(alpha = alpha, classlist = AAclass, Gamma_0 = 1.5, Gamma_1 = 0.05)
	tmp = recom_1
	tmp[recom_1 == -1] = 9
	NB_predict(tmp, theta, maxL = 9, maxR = 10)
}

