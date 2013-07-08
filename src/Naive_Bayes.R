###This is a file containing functions of the Naive Bayes model.
#Contents:
#		getFeatures:
#			Get the feature vectors of all peptides. 
#		Dirichlet_Parameter:
#			Computes the posterior parameter vector
#                       alpha = (alpha_1,alpha_2,...alpha_K) of 
#			Dirichlet distribution.
#		getTheta_MC:
#			A supporting function that computes the likelihood parameters of the 
#           Naive Bayesian model, assuming the likelihood has Dirichlet 
#    		distribution. The likelihood is simulated.
#		NB_predict:
#			Predict the outcome value of each peptide in newdata.
#
##===============================================================================




getFeatures <- function(data.org, AAclass, nL, nR)
{
#================================================================================
#Function: getFeatures
#
#--------------------------------------------------------------------------------
#Description:    
#    Get the feature vectors of all peptides. 
#
#--------------------------------------------------------------------------------
#Input arguments:
#    data.org
#             A matrix whose rows correspond to peptides. The matrix has n+2 
#             columns. The first n show whether each peptide works with n
#             different enzymes. The last two columns contain amino-acids at the
#             left of the serene('nterms') and right('cterms').
#    class
#             A matrix specifying the mapping between each amino-acids and the 
#             class it belongs to
#    nL
#             Number of amino acids at the left of the serene that will be 
#             features.
#    nR
#             Number of amino acids at the right of the serene that will be 
#             features.
#
#--------------------------------------------------------------------------------
#Return objects:
#    feature
#             A data frame whose rows correspond to peptides and columns to 
#             features. The last columns are the outcome values(whether each 
#             peptide works with enzyme 1, 2a and 2b).
#
#--------------------------------------------------------------------------------
    #nVal: number of values each feature can take
    nVal <- max(AAclass)
    #nOUTCOME: number of outcome values
    nOUTCOME <- dim(data.org)[2]-2
    feature <- matrix(NA, nrow=dim(data.org)[1], ncol=nL+nR+nOUTCOME)
    for (r in 1:dim(data.org)[1]) {
        sequence <- unlist(strsplit(data.org[r, 'nterm'],''))
        l.seq  <- length(sequence)
        l <- min(nL,l.seq)
        for (i in 1:l) {
		feature[r,nL+1-i] <- AAclass[1,sequence[l.seq+1-i]]
	  }
        sequence <- unlist(strsplit(data.org[r, 'cterm'],''))
        l.seq  <- length(sequence)
        c <- min(nR,l.seq)
        for (i in 1:c) {
		feature[r,nL+i] <- AAclass[1,sequence[i]]
	  }
        if( nOUTCOME != 0) {
            for (i in 1:nOUTCOME) {
                feature[r,nL+nR+i] <- data.org[r,i]
            }
        }
    }    
    #outcome.names : name of outcome values
    outcome.names <- c()
    if(nOUTCOME != 0) {
        outcome.names <- colnames(data.org)[1:nOUTCOME]
    }
    feature <- as.data.frame(feature)
    if(nOUTCOME != 0){
        colnames(feature) <- c(paste('L',nL:1,sep=""),
            paste('R',1:nR,sep=""), outcome.names)
    }
    else {
        colnames(feature) <- c(paste('L',nL:1,sep=""), 
            paste('R',1:nR,sep=""))
    }
    return( feature )
}




##==============================================================================
Dirichlet_Parameter <- function(train, nVal) 
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
#    nVal
#             Number of values each feature can take.
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




##==============================================================================
etTheta_MC <- function(train = NA, alpha = NA, nVal) 
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



##===============================================================================
NB_predict <- function(newdata, theta, prior.positive = 10**(-4))
{
#================================================================================
#Function: NB_predict
#
#--------------------------------------------------------------------------------
#Description:
#    Predict the outcome value of each peptide in newdata.
#
#--------------------------------------------------------------------------------
#Input arguments:
#    newdata
#            A data frame whose rows correspond to peptides and columns to 
#            features. THERE SHALL BE NO OUTCOME VALUES.
#    theta
#            Likelihood parameters. A list returned by the function getTheta.
#    prior.positive
#            Prior distribution Pr(Y=1) where Y is the outcome value. The default 
#            10^-4. 
#
#--------------------------------------------------------------------------------
#Return objects:
#    predict
#            A vector of posterior outcomes of each peptide.  
#--------------------------------------------------------------------------------
    theta_0 <- theta$theta_0
    theta_1 <- theta$theta_1
    nData <- dim(newdata)[1]
    #K: number of features
    K <- dim(newdata)[2] 
    predict <- rep(0, nData)
    for(n in 1:nData) {
        feature <- as.numeric(newdata[n,])
        likelihood_1 <- 1
        likelihood_0 <- 1
        for(i in 1:K) {
            if(!is.na(feature[i])){
                likelihood_1 <- theta_1[feature[i],i]*likelihood_1
                likelihood_0 <- theta_0[feature[i],i]*likelihood_0
            }
        }
        predict[n] <- likelihood_1*prior.positive/(likelihood_0*(1-prior.positive)
            + likelihood_1*prior.positive)
    }
    return(predict)
}