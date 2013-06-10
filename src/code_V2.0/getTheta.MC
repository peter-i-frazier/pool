getTheta.MC <- function(train, nVal) 
{
#===============================================================================
#Function: getTheta.MC
#
#-------------------------------------------------------------------------------
#Description:
#    A supporting function that computes the likelihood parameters of the Naive 
#    Bayesian model, assuming the likelihood has (posterior) Dirichlet distrib-
#    -ution. The likelihood is simulated.
# 
#-------------------------------------------------------------------------------
#Input arguments:
#    train
#             A list returned by the function getFeature. 
#
#-------------------------------------------------------------------------------
#Return objects:
#    theta
#             A list of two likelihood matrices: theta_0 and theta_1. Each matr-
#             -ix corresponds to a likelihood matrix. The rows of a likelihood 
#             matrix correspond to possible values a feature can take and colum-
#             -mns to features. Each entry correponds to a likelihood probability.
#
#-------------------------------------------------------------------------------
	
      # K: number of features
      K <- dim(train)[2] - 1
      outcome.value = train[,K+1]
      # devide train data by outcome value
      train_0 <- train[outcome.value == 0, 1:K]
      train_1 <- train[outcome.value == 1, 1:K]   
      
	R0 <- dim(train_0)[1]
      R1 <- dim(train_1)[1]
	theta_0 <- c()
      theta_1 <- c()

	for (col in 1:K) {
		count_0 <- rep(0,nVal)
		for (r in 1:R0) {
			count_0[train_0[r,col]] <- count_0[train_0[r,col]] + 1
			alpha_0 <- count_0 + rep(abs(col-(K+1)/2)**0.5,nVal)
		}
            gamma.rdn <- rgamma(nVal, shape = alpha_0)	
            theta_0 <- cbind(theta_0, gamma.rdn/sum(gamma.rdn))	
            count_1 <- rep(0,nVal)            
		for (r in 1:R1) {
			count_1[train_1[r,col]] <- count_1[train_1[r,col]] + 1
			alpha_1 <- count_1 + rep(abs(col-(K+1)/2)**0.5,nVal)
		}	
            theta_1 <- cbind(theta_1, alpha_1/sum(alpha_1))
            gamma.rdn <- rgamma(nVal, shape = alpha_1)	
            theta_1 <- cbind(theta_1, gamma.rdn/sum(gamma.rdn))	

	}
      theta = list("theta_0" = theta_0, "theta_1" = theta_1)
	return (theta)
}


