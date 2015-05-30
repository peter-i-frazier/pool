getTheta_MC <- function(train = NA, alpha = NA, classlist) 
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
#    classlist
#             A matrix specifying the mapping between each amino-acids and the 
#             class it belongs to.    
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
	nVal <- length(unique(as.numeric(classlist)))
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




####    This is the function: Dirichlet_Parameter required by getTheta_MC    ####
Dirichlet_Parameter <- function(train, classlist) 
{
#===============================================================================
#Function: Dirichlet_Parameter
#-------------------------------------------------------------------------------
#Description:
#    Computes the posterior parameter vector
#                       alpha = (alpha_1,alpha_2,...alpha_K) of 
#    Dirichlet distribution Dir(alpha,K) for each feature given training data. 
#    The prior parameter vector of a feature d positions away from the serene is
#    set as
#                       (d**0.5,...,d**0.5).
# 
#-------------------------------------------------------------------------------
#Input arguments:
#    train
#             A data matrix. Each row corresponds to a peptide and each column 
#             corresponds to a feature. The last column contain outcome values.
#    classlist
#             A matrix specifying the mapping between each amino-acids and the 
#             class it belongs to.
#
#-------------------------------------------------------------------------------
#Return objects:
#    alpha
#             A list of two data matrices specifying the posterior alpha 
#             parameter of Dirichlet distribution, with the outcome value being 0
#             and 1 respectively.
#             For each matrix, the columns correspond to features. Each column 
#             contains the alpha parameter of a feature.
#
#-------------------------------------------------------------------------------
    #nVal: number of values a feature can take
	nVal <- length(unique(as.numeric(classlist)))
	# K: number of features
    K <- dim(train)[2] - 1
    outcome.value = train[,K+1]
    # divide training data by outcome value
    train_0 <- train[outcome.value == 0, 1:K]
    train_1 <- train[outcome.value == 1, 1:K]   
	    
	R_1 <- dim(train_1)[1]
    alpha_1 <- c()
	for (col in 1:K) {
		count <- rep(0,nVal)
		for (r in 1:R_1) {
                  if (!is.na(train_1[r,col])) {
                      count[train_1[r,col]] <- count[train_1[r,col]] + 1
                  }
		}
		distance <- as.numeric(paste(unlist(strsplit(colnames(train_1)[col],''))[-1],collapse=""))
        alpha_1 <- cbind(alpha_1, count + rep(distance**0.5,nVal))
    }
    R_0<- dim(train_0)[1]
    alpha_0 <- c()
	for (col in 1:K) {
		count <- rep(0,nVal)
		for (r in 1:R_0) {
                  if (!is.na(train_0[r,col])) {
                      count[train_0[r,col]] <- count[train_0[r,col]] + 1
                  }
		}
		distance <- as.numeric(paste(unlist(strsplit(colnames(train_0)[col],''))[-1],collapse=""))
        alpha_0 <- cbind(alpha_0, count + rep(distance**0.5,nVal))
    }	
	alpha_0 <- as.data.frame(alpha_0)
    colnames(alpha_0) <- colnames(train)[1:K]
    alpha_1 <- as.data.frame(alpha_1)
    colnames(alpha_1) <- colnames(train)[1:K]	
    alpha <- list("alpha_0" = alpha_0, "alpha_1" = alpha_1)	
    return (alpha)
} 

