getTheta_MC <- function(train = NA, alpha = NA, nVal) 
{
#===============================================================================
#Function: getTheta.MC
#
#-------------------------------------------------------------------------------
#Description:
#    A supporting function that computes the likelihood parameters of the Naive 
#    Bayesian model, assuming the likelihood has Dirichlet 
#    distribution. The likelihood is simulated.
# 
#-------------------------------------------------------------------------------
#Input arguments:
#    train
#             A data matrix. Each row corresponds to a peptide and each column 
#             corresponds to a feature with the last column being the outcome
#             value.
#    alpha
#             A list of two data matrices specifying the posterior alpha 
#             parameter of Dirichlet distribution, with the outcome value being 
#             0 and 1 respectively.
#             For each matrix, the columns correspond to features. Each column 
#             contains the alpha parameter of a feature.       
#    train and alpha need not to be specified together. When both are specified, 
#    train will be ignored.
#
#    nVal
#             Number of values each feature can take.
#
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
	if(missing(alpha)){
	   alpha <- Dirichlet_Parameter(train, nVal)
	}
    alpha_0 <- alpha$alpha_0
	alpha_1 <- alpha$alpha_1
	K <- dim(alpha_1)[2]
	theta_0 <- c()
    theta_1 <- c()

	for (col in 1:K) {
        gamma.rdn <- rgamma(nVal, shape = alpha_0[,col])	
        theta_0 <- cbind(theta_0, gamma.rdn/sum(gamma.rdn))
            
        gamma.rdn <- rgamma(nVal, shape = alpha_1[,col])	
        theta_1 <- cbind(theta_1, gamma.rdn/sum(gamma.rdn))         
	}
    theta_0 <- as.data.frame(theta_0)
    colnames(theta_0) <- colnames(alpha_0)
    theta_1 <- as.data.frame(theta_1)
    colnames(theta_1) <- colnames(alpha_1)
    theta <- list("theta_0" = theta_0, "theta_1" = theta_1)
	return (theta)
}


