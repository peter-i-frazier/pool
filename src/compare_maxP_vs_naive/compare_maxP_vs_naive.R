rm(list=ls())

#=================================================================================
#Specify Paths and working directory
dataFile <- '../../data/newData.csv'
classFile <- '../../data/Reduced_AA_Alphabet.csv'
#=================================================================================
#import module
source('../Naive_Bayes/Naive_Bayes_util.R')
source('../Opt_Search/Opt_Search_util.R')

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
# rdn peptide library
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
print ('rdn peptide lib done')

# mutate library
Mut <- 1e4
mutlib <- c()
orig.hits <- trainData[trainData[,'AcpH']==1,1:38]
no.hits <- dim(orig.hits)[1]
for (i in 1:Mut) {
	# choose which peptide to mutate
	pep.toMut <- orig.hits[ceiling(no.hits*runif(1)),]
	# choose length of mutated peptide
	tL <- ceiling(runif(1)*(maxLen - minLen)) + minLen
	ratioL <- runif(1)*0.4 + 0.2
	L <- floor(tL * ratioL)
	R <- tL - L
	for (j in 1:length(pep.toMut)) {
		if (j<=(nF/2-L)) {
			pep.toMut[j] <- -1
		}
		if(j>=(nF/2+R+1)) {
			pep.toMut[j] <- -1
		}
	} 
	# locate start and end point
	for (j in 1:nF) {
		if (pep.toMut[j] != -1) {
			start <- j
			break
		}
	}
	for (j in start:nF) {
		if (pep.toMut[j] == -1) {
			end <- j-1
			break
		}
	}
	# choose no. of positions to mutate
	no.pos <- ceiling(runif(Nlib)*4)  
	# locate mutating position
	pos <- ceiling(runif(no.pos[n]) * (end-start+1)) + start - 1
	for (k in 1:length(pos)) {
		pep.toMut[pos[k]] <- ceiling(runif(1) * nAA)
	}
	mutlib <- rbind(mutlib, pep.toMut)
	# print (i)
}
mutlib <- as.matrix(mutlib)
print ('mutlib done')

# combined library
combo.lib <- rbind(peplib, mutlib)

maxL <- ceiling(maxLen*0.6)
maxR <- ceiling(maxLen*0.8)
#recom
print ('add in=10, includes mutate')
recom <- maxP_search3(X, Y, AAclass, Nrec=100, itr=400, combo.lib, maxL=maxL, maxR=maxR, Gamma_1 = Gamma_1, Gamma_0 = Gamma_0, add_ins = 10)

#Naive method
Itr <- 1000
no.rec <- 100
prob <- 0
alpha <- Dirichlet_Parameter(X, Y, AAclass, Gamma_0 = Gamma_0, Gamma_1 = Gamma_1)
for (i in 1:Itr) {
	theta <- getTheta_MC(alpha = alpha, classlist = AAclass, Gamma_0 = Gamma_0, Gamma_1 = Gamma_1)
	prob <- NB_predict(combo.lib, theta, maxL = maxL, maxR = maxR) + prob
}
prob <- prob/ Itr
dec_order <- order(prob, decreasing=T)
naive_rec <- combo.lib[dec_order[1:no.rec],]
naive_idx <- dec_order[1:no.rec]
naive_prob <- prob[dec_order[1:no.rec]]

# calculate prob for recom
recom_prob <- 0
for (i in 1:Itr) {
	theta <- getTheta_MC(alpha = alpha, classlist = AAclass, Gamma_0 = Gamma_0, Gamma_1 = Gamma_1)
	recom_prob <- NB_predict(recom$rec, theta, maxL = maxL, maxR = maxR) + recom_prob
}
recom_prob <- recom_prob/ Itr

