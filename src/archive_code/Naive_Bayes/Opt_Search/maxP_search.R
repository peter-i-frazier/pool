#================================================================================
#Function: maxP_search
#
#--------------------------------------------------------------------------------
#Description:
#    give a set of recommendations S that maximize P(at least one peptide in S is hit)
#
#--------------------------------------------------------------------------------
#Input arguments:
#	 data
#			 data.frame with 1st column response, and rest columns are features
#	 classlist
#             A matrix specifying the mapping between each amino-acids and the 
#             class it belongs to.
#    Nrec
# 			 No. of recommendations to be generated
# 	 itr
#			 Iteration needed for calculation prob being hit
#	 Nlib
# 			 No. of peptides considered for searching optimal
#	 maxL
#	 		 max no. of AAs to the left of Serine
#	 maxR
#	 		 max no. of AAs to the right of Serine
#	 rootpath
#	 		 path of src folder
#
#--------------------------------------------------------------------------------
#Return objects:
#    a list
#            result$rec:	matrix with each row representing a recommended peptide  
#--------------------------------------------------------------------------------
maxP_search <- function(data, classlist, Nrec, itr=1000, Nlib = 1e6, maxL, maxR, rootpath) {
	source(paste0(rootpath, '/Naive_Bayes/Naive_Bayes_util.R'))
	nAA <- length(unique(as.numeric(classlist)))
	# X.1 <- data[data[,'Y']==1,][,-1]
	# X.0 <- data[data[,'Y']==0,][,-1]
	Y <- data[,dim(data)[2]]
	class(Y)
	X <- data[,-dim(data)[2]]
	rec <- matrix(NA, nrow=Nrec, ncol=(maxL+maxR))
	for (n in 1:Nrec) {
		# sample random peptides
		L <- ceiling(runif(Nlib) * (maxL - 4)) + 4
		R <- ceiling(runif(Nlib) * (maxR - 4)) + 4
		rdn.peptides <- matrix(NA, nrow=Nlib, ncol=(maxL+maxR))
		for (j in 1:Nlib) {
			peptide <- ceiling(runif(L[j]+R[j]) * nAA)
			rdn.peptides[j, (maxL-L[j]+1):(maxL+R[j])] <- peptide
		}
		# train & predict
		prob <- 0
		for (i in 1:itr) {
			theta <- getTheta_MC(X, Y, classlist=classlist, Gamma = 0.3)
			prob <- NB_predict(rdn.peptides, theta) + prob
		}
		prob <- prob/ itr
		# find one peptide with highest prob being hit, is recommendation
		rec[n,] <- rdn.peptides[which(prob==max(prob)),]
		X <- rbind(X, rec[n,])
		Y <- c(Y, 0)
	}
	result <- list('rec'=rec)
	return (result)
}