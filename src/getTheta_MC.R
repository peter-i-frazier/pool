getTheta_MC <- function(train, nVal) 
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
#             A data matrix. Each row corresponds to a peptide and each column 
#             correpsonds to a feature with the last column being the outcome
#             value.
#    nVal
#             Number of values each feature can take.
#-------------------------------------------------------------------------------
#Return objects:
#    theta
#             A list of two likelihood matrices: theta_0 and theta_1. 
#             The rows of a likelihood matrix correspond to possible values a 
#             feature can take and columns to features. Each entry (x,y) 
#             corresponds to a likelihood probability that feature y takes value 
#             x.
#
#-------------------------------------------------------------------------------
	
      # K: number of features
      K <- dim(train)[2] - 1
      outcome.value = train[,K+1]
      # devide train data by outcome value
      train_0 <- train[outcome.value == 0, 1:K]
      train_1 <- train[outcome.value == 1, 1:K]   

      alpha_0 <- Dirichlet_Parameter(train_0, nVal)
      alpha_1 <- Dirichlet_Parameter(train_1, nVal)
     	theta_0 <- c()
      theta_1 <- c()

	for (col in 1:K) {
            gamma.rdn <- rgamma(nVal, shape = alpha_0[col])	
            theta_0 <- cbind(theta_0, gamma.rdn/sum(gamma.rdn))
            
            gamma.rdn <- rgamma(nVal, shape = alpha_1[col])	
            theta_1 <- cbind(theta_1, gamma.rdn/sum(gamma.rdn))
            
	}
      theta_0 <- as.data.frame(theta_0)
      colnames(theta_0) <- colnames(train)[1:K]	
      theta_1 <- as.data.frame(theta_1)
      colnames(theta_1) <- colnames(train)[1:K]		
      theta = list("theta_0" = theta_0, "theta_1" = theta_1)
	return (theta)
}


