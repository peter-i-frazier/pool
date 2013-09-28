##############					Naive_Bayes					##################
#	INPUT:
#		trainX
#		trainY
#		testData
#		prior.positive
#		classlist
#		predIter
#	OUTPUT:
#		predProb

Naive_Bayes <- function(trainX, trainY, testData, classlist, Gamma = 1, prior.positive = 1e-4, predIter = 1000)
{
	nF <- dim(trainX)[2];
	nAA <- length(unique(as.numeric(classlist))) 
	#theta_0 <- matrix(table(as.numeric(classlist))/length(classlist),nAA,nF)
	#colnames(theta_0) <- colnames(trainX)
	predict_mat <- c()
	for(i in 1:predIter) {
		theta <- getTheta_MC(trainX = trainX, trainY = trainY, classlist = AAclass, Gamma = Gamma)
		#theta$theta_0 <- theta_0
		if(is.vector(testData)){
			predict_mat <- c(predict_mat, NB_predict(testData, theta, prior.positive))
		}
		else{
			predict_mat <- rbind(predict_mat, NB_predict(testData, theta, prior.positive))
		}
	}
	predProb <- c()
	if(is.vector(testData)){
		predProb <- mean(predict_mat)
	}
	else{
		predProb <- colMeans(predProb)
	}
	return(predProb)
}




##############################################################################################
getFeatures <- function(data.org, classlist, nL, nR)
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
#    classlist
#             A matrix specifying the mapping between each amino-acids and the 
#             class it belongs to.
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
    nVal <- length(unique(as.numeric(classlist)))
    #nOUTCOME: number of outcome values
    nOUTCOME <- dim(data.org)[2]-2
    feature <- matrix(-1, nrow=dim(data.org)[1], ncol=nL+nR+nOUTCOME)
    for (r in 1:dim(data.org)[1]) {
        sequence <- unlist(strsplit(data.org[r, 'nterm'],''))
        l.seq  <- length(sequence)
        l <- min(nL,l.seq)
        for (i in 1:l) {
		feature[r,nL+1-i] <- classlist[1,sequence[l.seq+1-i]]
	  }
        sequence <- unlist(strsplit(data.org[r, 'cterm'],''))
        l.seq  <- length(sequence)
        c <- min(nR,l.seq)
        for (i in 1:c) {
		feature[r,nL+i] <- classlist[1,sequence[i]]
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

###############################################################################################
getTheta_MC <- function(trainX = NA, trainY = NA, alpha = NA, classlist, Gamma = 1) 
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
	   alpha <- Dirichlet_Parameter(trainX, trainY, classlist, Gamma)
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
Dirichlet_Parameter <- function(trainX, trainY, classlist, Gamma = 1) 
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
    K <- dim(trainX)[2]
    # divide training data by outcome value
    train_0 <- trainX[trainY == 0,]
    train_1 <- trainX[trainY == 1,]   
	
	R_1 <- dim(train_1)[1]
    alpha_1 <- c()
	for (col in 1:K) {
		count <- rep(0,nVal)
		for (r in 1:R_1) {
            if (train_1[r,col] != -1) {
                count[train_1[r,col]] <- count[train_1[r,col]] + 1
            }
		}
		distance <- as.numeric(paste(unlist(strsplit(colnames(train_1)[col],''))[-1],collapse=""))
        alpha_1 <- cbind(alpha_1, count + rep(distance**0.5*Gamma,nVal))
    }
    R_0<- dim(train_0)[1]
    alpha_0 <- c()
	for (col in 1:K) {
		count <- rep(0,nVal)
		for (r in 1:R_0) {
            if (train_0[r,col] != -1) {
                count[train_0[r,col]] <- count[train_0[r,col]] + 1
            }
		}
    	distance <- as.numeric(paste(unlist(strsplit(colnames(train_0)[col],''))[-1],collapse=""))
        alpha_0 <- cbind(alpha_0, count + rep(distance**0.5*Gamma,nVal))
    }	
	alpha_0 <- as.data.frame(alpha_0)
    colnames(alpha_0) <- colnames(trainX)
    alpha_1 <- as.data.frame(alpha_1)
    colnames(alpha_1) <- colnames(trainX)
    alpha <- list("alpha_0" = alpha_0, "alpha_1" = alpha_1)	
    return (alpha)
} 




#######################################################################################
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
	if(is.vector(newdata)){
		feature <- as.numeric(newdata)
        likelihood_1 <- 1
        likelihood_0 <- 1
		K <- length(newdata)
        for(i in 1:K) {
            if(feature[i]!=-1){
                likelihood_1 <- theta_1[feature[i],i]*likelihood_1
                likelihood_0 <- theta_0[feature[i],i]*likelihood_0
            }
        }
        predict <- likelihood_1*prior.positive/(likelihood_0*(1-prior.positive)
            + likelihood_1*prior.positive)
	}
	else{
		nData <- dim(newdata)[1]
		K <- dim(newdata)[2] 
		predict <- rep(0, nData)
		for(n in 1:nData) {
			feature <- as.numeric(newdata[n,])
			likelihood_1 <- 1
			likelihood_0 <- 1
			for(i in 1:K) {
				if(feature[i]!=-1){
					likelihood_1 <- theta_1[feature[i],i]*likelihood_1
					likelihood_0 <- theta_0[feature[i],i]*likelihood_0
				}
			}
			predict[n] <- likelihood_1*prior.positive/(likelihood_0*(1-prior.positive)
					+ likelihood_1*prior.positive)
		}
	}    
    return(predict)
}





#######################################################################################
NB_predict_par <- function(newdata, theta, prior.positive = 10**(-4))
{
#================================================================================
#Function: NB_predict_par
#
#--------------------------------------------------------------------------------
#Description:
#    Predict the outcome value of each peptide in newdata. This is a parallel 
#    version of NB_predict. When newdata is large(e.g, more than  peptides), 
#    this function should be preferred.
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
    predict <- foreach(n=1:nData, .init = c(), .combine = append) %dopar% {
        feature <- as.numeric(newdata[n,])
        likelihood_1 <- 1
        likelihood_0 <- 1
        for(i in 1:K) {
            if(feature[i]!=-1){
                likelihood_1 <- theta_1[feature[i],i]*likelihood_1
                likelihood_0 <- theta_0[feature[i],i]*likelihood_0
            }
        }
        likelihood_1*prior.positive/(likelihood_0*(1-prior.positive)
            + likelihood_1*prior.positive)
    }
    return(predict)
}




#################################################################################################
NB_peptideSim_par <- function(theta = NA, train = NA, outcome_name = NA, 
 classlist, nPep, isHit = 1, minL = 1, maxL = 3, minR = 3, maxR = 9 )
{
#================================================================================
#Function: NB_peptideSim_par
#
#--------------------------------------------------------------------------------
#Description:
#    This is the parallel version of NB_peptideSim. When nPep is large, this 
#    function should be preferred.
#         Simulate a set of peptides. either totally randomly or according to
#    theta. Theta can be specified or trained using training data. 
#    The length of simulated peptides:
#        Currently, the length to the left & the right of the serene are drawn 
#    from independent uniform distributions, where the left is uniform from 
#    length minL(default value is 1) to maxL(default value is 3), and the right 
#    is uniform from length minR(default value is 3) to maxR(default value is 9).
#        For further improvement, the distribution of the length of peptides 
#    could be a parameter of this function.
#
#--------------------------------------------------------------------------------
#Input arguments:
#    theta
#             A list of two likelihood matrices: theta_0 and theta_1. 
#             The rows of a likelihood matrix correspond to possible values a 
#             feature can take and columns to features. Each entry (x,y) 
#             corresponds to a likelihood probability that feature y takes value 
#             x. 
#    train
#             A data matrix. Each row corresponds to a peptide and each column 
#             corresponds to a feature with the last column being the outcome
#             value.
#    The user may either specify train or theta, or neither. If both are 
#    specified, TRAIN is ignored. If neither is specified, peptides will be
#    generated totally randomly.
#    outcome_name
#             Name of the outcome value. 
#    classlist
#             A matrix specifying the mapping between each amino-acids and the 
#             class it belongs to.
#    nPep
#             Number of peptides to be generated.
#    isHit
#             A binary parameter. "1" if the peptides are to be generated 
#             from theta parameters conditioning on positive outcome value, "0"
#             otherwise. The default value is 1.              
#
#--------------------------------------------------------------------------------
#Return objects:
#    peptides
#            A data matrix with 3 columns. Each line corresponds to a generated 
#            peptide. The first column is the outcome value. The last two 
#            correspond to nterms(amino acids at the left of the serene)and 
#            cterms(amino-acids at the right of the serene).
#--------------------------------------------------------------------------------
nVal <- length(unique(as.numeric(classlist)))
#Case 1 : Generate peptides totally randomly
if(missing(train) && missing(theta)) {
	peptides <- foreach(icount(nPep), .combine = rbind, .init = c()) %dopar% {
		#Draw left & right lengths at random from uniform distribution
		sL <- minL-1+ceiling((maxL-minL+1)*runif(1))
		sR <- minR-1+ceiling((maxR-minR+1)*runif(1))
		l_peptide <- c()
		r_peptide <- c()
		#Iterate through positions
		for(i in sL:1) {			
			s_class <- ceiling(8*runif(1)) 
			s_class_list <- colnames(classlist)[classlist == s_class]
			s_amino_acid <- s_class_list[ceiling(length(s_class_list)*runif(1))]
			l_peptide <- paste(l_peptide, s_amino_acid, sep="")
		}
		for(i in 1:sR) {		
			s_class <- ceiling(8*runif(1)) 
			s_class_list <- colnames(classlist)[classlist == s_class]
			s_amino_acid <- s_class_list[ceiling(length(s_class_list)*runif(1))]
			r_peptide <- paste(r_peptide, s_amino_acid, sep="")
		}   
    c(l_peptide,r_peptide)
	}
}
#Case 2: Generate peptides according to theta or training data(if theta is not specified)
else {
    alpha <- NA
    if(missing(theta)) {
	alpha <- Dirichlet_Parameter(train,classlist)
	}
	peptides <- foreach(icount(nPep), .init = c(), .combine = rbind, 
	    .multicombine = TRUE, .export = c('rdist','getTheta_MC')) %dopar% {
		#Draw left & right lengths at random from uniform distribution
		sL <- minL-1+ceiling((maxL-minL+1)*runif(1))
		sR <- minR-1+ceiling((maxR-minR+1)*runif(1))
	    if(!is.na(alpha)){
			#Simulate theta according to alpha parameter
			theta <- getTheta_MC(alpha = alpha, classlist = classlist)	
	    }
        if(isHit) {
	        param <- theta$theta_1 }
	    else {
			param <- theta$theta_0 }		
		l_peptide <- c()
		r_peptide <- c()
		#Iterate through positions
		for(i in sL:1) {
			s_class <-  rdist(c(1:nVal), param[,paste('L',i,sep="")])
		    s_class_list <- colnames(classlist)[classlist == s_class]
			s_amino_acid <- s_class_list[ceiling(length(s_class_list)*runif(1))]
			l_peptide <- paste(l_peptide, s_amino_acid, sep="")
		}
		for(i in 1:sR) {		
			s_class <-  rdist(c(1:nVal), param[,paste('R',i,sep="")])
			s_class_list <- colnames(classlist)[classlist == s_class]
			s_amino_acid <- s_class_list[ceiling(length(s_class_list)*runif(1))]
			r_peptide <- paste(r_peptide, s_amino_acid, sep="")
		}   
		c(l_peptide,r_peptide)
	}
}
result <- data.frame(matrix(-1, nrow = nPep, ncol = 3))
result[,1] <- as.numeric(rep(isHit, nPep))
result[,-1] <- peptides
colnames(result) <- c(outcome_name,'nterm','cterm')
return(result)
}



##########    This is the function: rdist required by NB_peptideSim    #########
rdist <- function(value, p) 
{
#===============================================================================
#Function: rdist
#
#-------------------------------------------------------------------------------
#Description:
#   Generate a random number that has a discrete distribution.
# 
#-------------------------------------------------------------------------------
##Input arguments:
#   value
#             A vector contain the possiblly taken values of the random number.
#   p
#             A vector contain the probability the random number taking each 
#             value. p must have a 1-1 correspondence with value.
#-------------------------------------------------------------------------------
#Return objects:
#   rdn
#             Generated random number.             
#
#-------------------------------------------------------------------------------
    l <- length(p)
    sum_p <- rep(0, l)
    for(i in 1:l) {
        sum_p[i] <- sum(p[1:i])
    }
    u_rdn <- runif(1)
    for(i in 1:l) {
        if(u_rdn < sum_p[i]){
            rdn <- value[i]
            break
        }
    }
    return(rdn)
}