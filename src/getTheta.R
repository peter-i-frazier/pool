getTheta <- function(train, nVal) 
{
#===============================================================================
#Function: getTheta
#-------------------------------------------------------------------------------
#Description:
#    A supporting function that computes the likelihood parameters of the Naive 
#    Bayesian model, assuming the likelihood has (posterior) Dirichlet 
#    distribution. The likelihood is the mean of the Dirichlet distribution.
# 
#-------------------------------------------------------------------------------
#Input arguments:
#    train
#             A data matrix. Each row corresponds to a peptide and each column 
#             corresponds to a feature with the last column being the outcome
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
    # divide train data by outcome value
    train_0 <- train[outcome.value == 0, 1:K]
    train_1 <- train[outcome.value == 1, 1:K]   
      
	R0 <- dim(train_0)[1]
    R1 <- dim(train_1)[1]
	theta_0 <- c()
    theta_1 <- c()

	for (col in 1:K) {
		count_0 <- rep(0,nVal)
		for (r in 1:R0) {
            if (!is.na(train_0[r,col])) {
                count_0[train_0[r,col]] <- count_0[train_0[r,col]] + 1
            }
		}
        alpha_0 <- count_0 + rep(abs(col-(K+1)/2)**0.5,nVal)	
        theta_0 <- cbind(theta_0, alpha_0/sum(alpha_0))
        count_1 <- rep(0,nVal)
		for (r in 1:R1) {
            if(!is.na(train_1[r,col])) {
		        count_1[train_1[r,col]] <- count_1[train_1[r,col]] + 1
            }
		}	
        alpha_1 <- count_1 + rep(abs(col-(K+1)/2)**0.5,nVal)
        theta_1 <- cbind(theta_1, alpha_1/sum(alpha_1))
	}
    theta = list(theta_0 = theta_0, theta_1 = theta_1)
	return (theta)
}