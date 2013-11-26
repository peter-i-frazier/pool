rm(list=ls())
source('../Naive_Bayes/Naive_Bayes_util.R')
#=================================================================================
# Utility functions
get_length <- function(peptide) {
	len <- 1
	for (i in 1:length(peptide)) {
		if (peptide[i] != -1){
			len <- len + 1
		}
	}
	return (len)
}
prob_shortest_hit <- function(peptides, prob, b) {
# Calculate P(shortest hit <= b)
	idx <- c()
	no <- 0
	if (is.vector(peptides)) {
		if (get_length(peptides)<=b){
			result <- prob
		} else {
			result <- 0
		}
	} else {
		for (n in 1:dim(peptides)[1]) {
			if (get_length(peptides[n,])<=b) {
				idx <- c(idx, n)
				no <- no +1
			}
		}
		if (no == 0) {
			result <- 0
		} else {
			prod <- 1
			for (i in 1:length(idx)) {
				prod <- prod * (1.0-prob[idx[i]])
			}
			result <- 1.0-prod
		}
	}
	return (result)
}

#=================================================================================
load('recom_data#2Benchmark.RData')
classFile <- '../../data/Reduced_AA_Alphabet.csv'
dataFileB2 <- '../../data/newData#2.csv'
data_orgB2 <- data.frame(read.csv(dataFileB2, header = T, as.is = T, sep = ","))
AAclass <- read.csv(classFile, header=T, as.is = T, sep=",")
#Set parameters
nL <- 19
nR <- 19
trainDataB2 <- getFeatures(data_orgB2,AAclass,nL,nR)
nAA <- max(AAclass[1,])
outcome <- 'AcpH'
XB2 <- as.matrix(trainDataB2[,c(1:(nL+nR))])
YB2 <- trainDataB2[,outcome]
Gamma_0 <- 1000
Gamma_1 <- 0.05
N <- 100
b <- 12
Itr <- 1000
# calculate posterior of theta
alpha <- Dirichlet_Parameter(XB2, YB2, AAclass, Gamma_0 = Gamma_0, Gamma_1 = Gamma_1)
# calculate prob_shortest_hit for recom
P.recom <- rep(0,N)
for (itr in 1:Itr) {
	theta <- getTheta_MC(alpha = alpha, classlist = AAclass, Gamma_0 = Gamma_0, Gamma_1 = Gamma_1)
	prob <- NB_predict(recom$rec, theta, maxL = nL, maxR = nR)
	P_short_hit <- rep(0,N)
	for (i in 1:N) {
		P_short_hit[i] <- prob_shortest_hit(recom$rec[1:i,], prob[1:i], b)
	}
	P.recom <- P.recom + P_short_hit
}
P.recom <- P.recom/Itr
print ('P.recom calculated')
# calculate prob_shortest_hit for naive
P.naive <- rep(0,N)
for (itr in 1:Itr) {
	theta <- getTheta_MC(alpha = alpha, classlist = AAclass, Gamma_0 = Gamma_0, Gamma_1 = Gamma_1)
	prob <- NB_predict(naive_rec, theta, maxL = nL, maxR = nR)
	P_short_hit <- rep(0,N)
	for (i in 1:N) {
		P_short_hit[i] <- prob_shortest_hit(naive_rec[1:i,], prob[1:i], b)
	}
	P.naive <- P.naive + P_short_hit
}
P.naive <- P.naive/Itr
print ('P.naive calculated')
# calculate prob_shortest_hit for mutate
P.mutate <- rep(0,N)
for (itr in 1:Itr) {
	theta <- getTheta_MC(alpha = alpha, classlist = AAclass, Gamma_0 = Gamma_0, Gamma_1 = Gamma_1)
	prob <- NB_predict(mutlib, theta, maxL = nL, maxR = nR)
	P_short_hit <- rep(0,N)
	for (i in 1:N) {
		for (j in 1:100) {
			P_short_hit[i] <- P_short_hit[i] + prob_shortest_hit(mutlib[(100*(j-1)+1):(100*(j-1)+i),], prob_recom[(100*(j-1)+1):(100*(j-1)+i)], b)
		}
	}
	P.mutate <- P.mutate + P_short_hit/100
}
P.mutate <- P.mutate/Itr
print ('P.mutate calculated')

# save P to csv for plot in MATLAB
PP <- rbind(P.recom, P.naive, P.mutate)
write.csv(PP, 'PP_use_newData#2.csv')



