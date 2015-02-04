# Given training data(trainX and trainY), train a model(learn alpha parameters) and predict the posterior Pr(Y=1|train data,alpha) for all points in testData, where Y is the label. The prediction is done by Monte Carlo simulation.	 
# 
# :param trainX: A matrix of training samples. Each row is the feature vector of a peptide, namely it is the first (nL + nR) columns of 'feature' in the function (see the note of function GetFeatures). 
# :param trainY: A vector of training labels. The length of trainY must equal the number of rows of trainX.
# :param testData: Data points whose labels are to be predicted. Arranged in the same format as trainX.
# :param classlist: The vector specifying the class each amino acid belongs to.    
# :param S.Pos: The position of serene in each feature vector. S.Pos is determined by the model, i.e., the theta matrix, if, for the model, there're at most nL amino-acids to the left of serene(i.e., there're nL theta parameters to the left of serene), then S.Pos = nL.	, maxR
# :param maxL: maxL is the number of peptides to the left of the serene of the peptide in testData with most number of peptides to the left of serene.
# :param maxR: same as maxL, except it's to the right of Serene
# :param gamma_0, gamma_1: The gamma parameters for computing posterior alpha parameters. When alpha is not specified, Gamma_0 and Gamma_1 must be specified.
# :param prior.positive: Prior probability  Pr(Y=1). Y is the label, default 10^-4.
# :param predIter: Number of iterations of simulation. In each iteration, a set of new theta parameter is simulated, and posterior Pr(Y=1|theta) is computed. 
# :return predProb: A vector of posterior probability of Pr(Y=1|train data) of each peptide in testData.  

Naive_Bayes <- function(trainX, trainY, testData, classlist, S.Pos, maxL, maxR, Gamma_0, Gamma_1, prior.positive, predIter)
{
	nF <- dim(trainX)[2]
	nAA <- length(unique(as.numeric(classlist))) 
	alpha <- Dirichlet_Parameter(trainX, trainY, classlist, Gamma_0 = Gamma_0, Gamma_1 = Gamma_1)
	predict_mat <- c()
	
	for(i in 1:predIter) {
		theta <- getTheta_MC(alpha = alpha, classlist = classlist, Gamma_0 = Gamma_0, Gamma_1 = Gamma_1)
		if(is.vector(testData)){
			predict_mat <- c(predict_mat, NB_predict(testData, theta, S.Pos, maxL, maxR, prior.positive))
		}
		else{
			predict_mat <- rbind(predict_mat, NB_predict(testData, theta, S.Pos, maxL, maxR, prior.positive))
		}
	}
	predProb <- c()
	if(is.vector(testData)){
		predProb <- mean(predict_mat)
	}
	else{
		predProb <- colMeans(predict_mat)
	}
	return(predProb)
}

NB_theta_matrices <- function(trainX, trainY, classlist, S.Pos, maxL, maxR, Gamma_0, Gamma_1, prior.positive)
{
	nF <- dim(trainX)[2]
	nAA <- length(unique(as.numeric(classlist))) 
	alpha <- Dirichlet_Parameter(trainX, trainY, classlist, Gamma_0 = Gamma_0, Gamma_1 = Gamma_1)
	
    theta <- getTheta(alpha = alpha, classlist = classlist, Gamma_0 = Gamma_0, Gamma_1 = Gamma_1)
	return(theta)
}

Naive_Bayes_new <- function(X_prefer, Y_prefer, X_unprefer, Y_unprefer, X_unlabel, Y_unlabel, testData, classlist, S.Pos, maxL, maxR, gamma_0_prefer, gamma_1_prefer, prior_prefer, gamma_0_unprefer, gamma_1_unprefer, prior_unprefer, gamma_0_unlabel, gamma_1_unlabel, prior_unlabel, itr) {
    alpha_prefer <- Dirichlet_Parameter(X_prefer, Y_prefer, classlist, gamma_0_prefer, gamma_1_prefer)
    alpha_unprefer <- Dirichlet_Parameter(X_unprefer, Y_unprefer, classlist, gamma_0_unprefer, gamma_1_unprefer)
    alpha_unlabel <- Dirichlet_Parameter(X_unlabel, Y_unlabel, classlist, gamma_0_unlabel, gamma_1_unlabel)

    predict_mat <- c()
    for (i in 1:itr) {
        theta_prefer <- getTheta_MC(alpha=alpha_prefer, classlist=classlist)
        theta_unprefer <- getTheta_MC(alpha=alpha_unprefer, classlist=classlist)
        theta_unlabel <- getTheta_MC(alpha=alpha_unlabel , classlist=classlist)
        predict_mat <- rbind(predict_mat, new_NB_predict(testData, theta_prefer, theta_unprefer, theta_unlabel, S.Pos, maxL, maxR, prior_prefer, prior_unprefer, prior_unlabel))
    }
    return (colMeans(predict_mat))
}


#	A supporting function that estimates the theta parameters of the Naive Bayesian model through Monte Carlo simulation, assuming theta has a Dirichlet distribution with 
#	hyperprior alpha. If alpha(posterior) is specified, theta will be simulated according to Dirichlet(alpha); otherwise trainX and trainY must be specified and posterior
#      	alpha will be estimated accoring to trainX and trainY by calling function Dirichlet_Parameter.
#===================================== 
#Input arguments:
#	trainX
#		A matrix of training samples. Each row is the feature vector of a peptide, namely it is the first (nL + nR) columns of feature(see the note of function
#       	getFeatures). 
#	trainY
#		A vector of training labels. The length of trainY must equal the number of rows of trainX.
#	alpha
#		A list of two data matrices specifying the alpha parameter of Dirichlet distribution:alpha_0 and alpha_1. Each is a nC*(nL+nR) matrix, where nC is the number of
#		classes, nL and nR has the same meaning in getFeatures. For each matrix, each column is the alpha parameter of the corresponding position.
#		E.g. alpha_0 has the followng format if nC = 5, nR = 5 and nC = 2:
#			1		2		...	10
#		1	alpha_{0,1}(1)	alpha_{0,2}(1)	...	alpha_{0,10}(1)
#		2	alpha_{0,1}(2)	alpha_{0,2}(2)	...	alpha_{0,10}(2)		
#	classlist
#		The vector specifying the class each amino acid belongs to.    
#	Gamma_0,Gamma_1
#		The gamma parameters for computing prior alpha parameters. When alpha is not specified, Gamma_0 and Gamma_1 must be specified.
#=====================================
#Return objects:
#	theta
#	A list of two matrices: theta_0 and theta_1, specifying the theta parameter. Each matrix has the same size as alpha_0/alpha_1 and is arranged in the same way as them. 
#=====================================       
getTheta_MC <- function(trainX = NA, trainY = NA, alpha = NA, classlist, Gamma_0 = NA, Gamma_1 = NA) 
{
	nVal <- length(unique(as.numeric(classlist)))
	if(missing(alpha)){
	   alpha <- Dirichlet_Parameter(trainX, trainY, classlist, Gamma_0, Gamma_1)
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

#	Similar as getTheta_MC. But theta is estimated as the mean of the Dirichlet distribution it follows. I.e., for k in {1..K}, where K is the number of classes,
#			thete_k = alpha_k/\Sum_{k=1}^K(alpha_k).
#===================================== 
#Input arguments:
#	see description of getTheta_MC
#=====================================
#Return objects:
#	see description of getTheta_MC
#=====================================
getTheta <- function(trainX = NA, trainY = NA, alpha = NA, classlist, Gamma_0, Gamma_1) 
{
       	nVal <- length(unique(as.numeric(classlist)))
	if(missing(alpha)){
	   alpha <- Dirichlet_Parameter(trainX, trainY, classlist, Gamma_0, Gamma_1)
	}
       	alpha_0 <- alpha$alpha_0
	alpha_1 <- alpha$alpha_1
	K <- dim(alpha_1)[2]
	theta_0 <- c()
	theta_1 <- c()
  
	for (k in 1:K) {	
        	theta_0 <- cbind(theta_0, alpha_0[,k]/sum(alpha_0[,k]))
        	theta_1 <- cbind(theta_1, alpha_1[,k]/sum(alpha_1[,k]))         
	}
	theta_0 <- as.data.frame(theta_0)
    	colnames(theta_0) <- colnames(alpha_0)
    	theta_1 <- as.data.frame(theta_1)
    	colnames(theta_1) <- colnames(alpha_1)
    	theta = list(theta_0 = theta_0, theta_1 = theta_1)
	return (theta)
}

#	This function estimates the posterior alpha parameter of the Dirichlet distribution theta follows. For the jth position to the left/right of serene of the feature vector
#	the prior alpha_1 is defined as(assuming K classes):
#		(Gamma_1*(j**0.5), Gamma_2*(j**0.5),...Gamma_2*(j**0.5))
#	and the prior alpha_0 is defined as(assuming K classes):
#		(Gamma_0*n_1/n, Gamma_0*n_2/n, ..., Gamma_0*n_K/n).
#	where for k in {1..K} n_k is the number of unique amino-acids in class k in the classlist and n is the total number of amino-acids in the classlist.
#===================================== 
#Input arguments:
#	trainX
#		see description of trainX in function getTheta_MC. 
#	trainY
#		see dexcription of trainY in function getTheta_MC.
#	classlist
#		see description of classlist in function getTheta_MC.
#	Gamma_0, Gamma_1
#		The gamma parameters for computing prior alpha parameters.
#=====================================
#Return objects:
#	alpha
#		see description of alpha in function getTheta_MC.   
#=====================================
Dirichlet_Parameter <- function(trainX, trainY, classlist, Gamma_0, Gamma_1) 
{
	#nVal: number of values a feature can take
	nVal <- length(unique(as.numeric(classlist)))
	# K: number of features
	K <- dim(trainX)[2]
	# divide training data by outcome value
	train_0 <- trainX[trainY == 0,]
	train_1 <- trainX[trainY == 1,]   
	
	#R_1: number of positive training samples
	R_1 <- dim(train_1)[1]
	alpha_1 <- c()
	for (col in 1:K) {
		count <- rep(0,nVal)
		for (r in 1:R_1) {
            if (train_1[r,col]!= -1) {
                count[train_1[r,col]] <- count[train_1[r,col]] + 1
            }
		}
		distance <- as.numeric(paste(unlist(strsplit(colnames(train_1)[col],''))[-1],collapse=""))
        alpha_1 <- cbind(alpha_1, count + rep(distance**0.5*Gamma_1,nVal))
    	}	
	
	#R_0: number of negative training samples
	R_0<- dim(train_0)[1]
    alpha_0 <- c()
	for (col in 1:K) {
		count <- rep(0,nVal)
		for (r in 1:R_0) {
            if (train_0[r,col]!= -1) {
                count[train_0[r,col]] <- count[train_0[r,col]] + 1
            }
		}
    		#distance <- as.numeric(paste(unlist(strsplit(colnames(train_0)[col],''))[-1],collapse=""))
        	alpha_0 <- cbind(alpha_0, count + table(as.numeric(classlist))/length(classlist)*Gamma_0)
    	}	
	alpha_0 <- as.data.frame(alpha_0)
    	colnames(alpha_0) <- colnames(trainX)
    	alpha_1 <- as.data.frame(alpha_1)
    	colnames(alpha_1) <- colnames(trainX)
    	alpha <- list("alpha_0" = alpha_0, "alpha_1" = alpha_1)	
    	return (alpha)
} 

#    	Predict label of each peptide in newdata. To speed up, when newdata contains more than 1 data point, we use the following way to compute predict:
#	For a data point x, assuming feature vector has length J, 
#		P(Y(x)=1|trainData, theta) = \Prod_{j=1}^J\theta_{1,j}(x_j)*Pr(Y=1) / (\Prod_{j=1}^J\theta_{1,j}(x_j)*Pr(Y=1) + \Prod_{j=1}^J\theta_{0,j}(x_j)*Pr(Y=0))
#					   = 1 / (1 + \Prod_{j=1}^J\eta_j(x_j)*(Pr(Y=0)/Pr(Y=1)))
#	\eta_j(x_j) = \theta_{0,j}(x_j) / \theta_{1,j}(x_j) is the ratio of theta_0 and theta_1.
#	A matrix ratio = theta_0/theta_1 is defined to store the eta matrix.	
#       In the case x_j = -1, \eta_j(-1) is undefined. To solve this problem, first append one line with all entries being 1 to the ratio matrix, now suppose the ratio matrix 
#	has m lines, then convert all newdata[n,j] to m for all n and j such that newdata[n,j] = -1. In this way when we multiplicate the ratio in all positions the result will
#	not be influenced by '-1' entries in newdata.	
#=====================================
#Input arguments:
#    	newdata
#            	Data points whose labels are to be predicted. Arranged in the same format as trainX in getTheta_MC.
#    	theta
#       	see the description of theta in the function.
#    	prior.positive
#            	Prior probability  Pr(Y=1). Y is the label, default 10^-4.
#	S.Pos
#		The position of serene in each feature vector. S.Pos is determined by the model, i.e., the theta matrix, if, for the model, there're at most nL amino-acids to 
#		the left of serene(i.e., there're nL theta parameters to the left of serene), then S.Pos = nL.	
#	maxL, maxR
#		maxL is the number of peptides to the left of the serene of the peptide in newdata with most number of peptides to the left of serene. Similar definition of
#       	maxR.
#=====================================
#Return objects:
#    	predict
#            A vector of posterior probability of Pr(Y=1|theta) of each peptide in newdata.  
#=====================================
NB_predict <- function(newdata, theta, S.Pos, maxL, maxR, prior.positive)
{
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
        	predict <- likelihood_1*prior.positive/(likelihood_0*(1-prior.positive) + likelihood_1*prior.positive)
	}
	else{
		nData <- dim(newdata)[1]
		nF <- dim(newdata)[2]
		newdata[newdata == -1] <- dim(theta_0)[1] + 1
		ratio <- theta_0/theta_1
		ratio <- rbind(ratio, rep(1, nF))
		colprod <- 1
		for(i in (S.Pos-maxL+1):(S.Pos+maxR)){ 
			colprod <- colprod*ratio[newdata[,i], i]
		}
		coeff <- (1-prior.positive)/prior.positive
		predict <- 1/(1+coeff*colprod)
	}    
    	return(predict)
}

new_NB_predict <- function(newdata, label_prefer_theta, label_unprefer_theta, unlabel_theta, S.Pos, maxL, maxR, label_prefer_prior, label_unprefer_prior, unlabel_prior) {
    return (NB_predict(newdata, label_prefer_theta, S.Pos, maxL, maxR, label_prefer_prior) * (1-NB_predict(newdata, label_unprefer_theta, S.Pos, maxL, maxR, label_unprefer_prior)) * NB_predict(newdata, unlabel_theta, S.Pos, maxL, maxR, unlabel_prior))
}

AUC <- function(x, y)
{
	n <- length(y)
	U <- 0
	L <- 0
	for(i in 1:(n-1)) {
		U <- U + y[i+1]*(x[i+1]-x[i])
		L <- L + y[i]*(x[i+1]-x[i])
	}
	return((U+L)/2)
}

# generate a single peptide
old_opt_gen_peptide_lib <- function(nF, maxL, maxR, minL, minR, S.Pos, alpha, classlist) {
    peptide <- rep(-1, nF)
    nL <- runif(1, min = minL, max = maxL)
    nR <- runif(1, min = minR, max = maxR)
    theta <- getTheta_MC(alpha=alpha, classlist=classlist)

    ratio <- as.matrix(theta$theta_1) / as.matrix(theta$theta_0)
    for (l in 1:nL) {
        peptide[S.Pos-l+1] <- which.max(ratio[,S.Pos-l+1])
    }
    for (r in 1:nR) {
        peptide[S.Pos+r] <- which.max(ratio[,S.Pos+r])
    }
    return (peptide)
}

# generate a single peptide
new_opt_gen_peptide <- function(nF, maxL, maxR, minL, minR, S.Pos, alpha_label_like, alpha_label_unlike, alpha_unlabel, classlist) {
    peptide <- rep(-1, nF)
    nL <- runif(1, min = minL, max = maxL)
    nR <- runif(1, min = minR, max = maxR)
    theta_label_like <- getTheta_MC(alpha=alpha_label_like, classlist=classlist)
    theta_label_unlike <- getTheta_MC(alpha=alpha_label_unlike, classlist=classlist)
    theta_unlabel <- getTheta_MC(alpha=alpha_unlabel, classlist=classlist)

    ratio <- as.matrix(theta_label_like$theta_1) / as.matrix(theta_label_like$theta_0) * as.matrix(theta_label_unlike$theta_0) / as.matrix(theta_label_unlike$theta_1) * as.matrix(theta_unlabel$theta_1) / as.matrix(theta_unlabel$theta_0)
    for (l in 1:nL) {
        peptide[S.Pos-l+1] <- which.max(ratio[,S.Pos-l+1])
    }
    for (r in 1:nR) {
        peptide[S.Pos+r] <- which.max(ratio[,S.Pos+r])
    }
    return (peptide)
}

# select the peptide to be added to recommendation set is not added before
select_new_recom <- function(trainX, peptides.library, prob) {
    unique_training_data <- unique(trainX)
    num_unique <- dim(unique_training_data)[1]
    while(length(prob)>0) {
        best_index <- which(prob==max(prob))
        add_pep <- peptides.library[best_index,]
        if((dim(unique(rbind(unique_training_data, add_pep)))[1] - num_unique) == 1){
            return (list('peptide'=add_pep, 'prob'=prob[best_index]))
        }
        else {
            peptides.library <- peptides.library[-best_index,]
            prob <- prob[-best_index]
        }
    }		
    return (rep(-10,1000))
}

get_ratio_old_method <- function(matrix_x, vector_y, classlist, gamma_0, gamma_1) {
    alpha <- Dirichlet_Parameter(matrix_x, vector_y, classlist, Gamma_0 = gamma_0, Gamma_1 = gamma_1) 
    theta <- getTheta(alpha = alpha, classlist = classlist, Gamma_0 = gamma_0, Gamma_1 = gamma_1)
    return (as.matrix(theta$theta_1) / as.matrix(theta$theta_0))
}

get_ratio_new_method <- function(classlist, matrix_x_label_prefer, vector_y_label_prefer, gamma_0_label_prefer, gamma_1_label_prefer, matrix_x_label_unprefer, vector_y_label_unprefer, gamma_0_label_unprefer, gamma_1_label_unprefer, matrix_x_unlabel, vector_y_unlabel, gamma_0_unlabel, gamma_1_unlabel) {
    alpha_label_prefer <- Dirichlet_Parameter(matrix_x_label_prefer, vector_y_label_prefer, classlist, Gamma_0 = gamma_0_label_prefer, Gamma_1 = gamma_1_label_prefer) 
    theta_label_prefer <- getTheta(alpha = alpha_label_prefer, classlist = classlist, Gamma_0 = gamma_0_label_prefer, Gamma_1 = gamma_1_label_prefer)

    alpha_label_unprefer <- Dirichlet_Parameter(matrix_x_label_unprefer, vector_y_label_unprefer, classlist, Gamma_0 = gamma_0_label_unprefer, Gamma_1 = gamma_1_label_unprefer) 
    theta_label_unprefer <- getTheta(alpha = alpha_label_unprefer, classlist = classlist, Gamma_0 = gamma_0_label_unprefer, Gamma_1 = gamma_1_label_unprefer)

    alpha_unlabel <- Dirichlet_Parameter(matrix_x_unlabel, vector_y_unlabel, classlist, Gamma_0 = gamma_0_unlabel, Gamma_1 = gamma_1_unlabel) 
    theta_unlabel <- getTheta(alpha = alpha_unlabel, classlist = classlist, Gamma_0 = gamma_0_unlabel, Gamma_1 = gamma_1_unlabel)
    return (as.matrix(theta_label_prefer$theta_1) / as.matrix(theta_label_prefer$theta_0) * as.matrix(theta_label_unprefer$theta_0) / as.matrix(theta_label_unprefer$theta_1) * as.matrix(theta_unlabel$theta_1) / as.matrix(theta_unlabel$theta_0))
}


isNew <- function(collection, to_add) {
# Test if to_add is already existed in collection. We only compare to_add's non -1 entries so
# if collection has "1 2 3 4 5 6 7" and to_add is "-1 -1 3 4 5 -1 -1" then we still return false 
# Inputs:
#   collection: matrix with more than 1 row
#   to_add: vector
# Output:
#   bool
    to_add_reconstruct <- c()
    collection_reconstruct <- c()
    for (j in 1:length(to_add)) {
        if (to_add[j] != -1) {
            to_add_reconstruct <- c(to_add_reconstruct, to_add[j])
            collection_reconstruct <- cbind(collection_reconstruct, collection[,j])
        }
    }
    return ((dim(unique(rbind(collection_reconstruct, to_add_reconstruct)))[1] - dim(unique(collection_reconstruct))[1]) == 1)
}

ROC_xy_output <- function(prob, Y) {
	FPR <- rep(-1, length(prob))
	TPR <- rep(-1, length(prob))
	thresholds <- sort(prob)
	for( i in 1:length(prob) ) {
		threshold <- thresholds[i]
		label <- rep(0, length(prob))
		for( j in 1:length(prob) ){
			if(prob[j] >= threshold) {
				label[j] <- 1 }
		}
		FPR[i] <- sum((Y==0)&(label==1))/sum(Y==0)
		TPR[i] <- sum((Y==1)&(label==1))/sum(Y==1)
	}
    return (list(x=FPR, y=TPR))
}

#	This function takes a matrix as input, each row is the feature vector of a peptide. It generates the amino-acid at each position uniformly according to classlist, and
#      	write generated peptides into a csv file and a text file named by out_file. 
#=====================================
#Input Arguments
#	newPep
#		The matrix of peptides to write, of which each row is the feature vector of the peptide.
#	S.Pos
#		The position of serene in each feature vector. S.Pos is determined by the model, i.e., the theta matrix, if, for the model, there're at most nL amino-acids to 
#		the left of serene(i.e., there're nL theta parameters to the left of serene), then S.Pos = nL.	
#	pred_prob
#		The vector of the predicted probability of being hit of all peptides in newPep.
#	classlist
#		The vector specifying the class each amino acid belongs to.    
#	out_file:
#		The name(including the full path) of output file. A csv file and a plane text file with this name will be generated to record the recommendations generated.
#=====================================
##	No objects are returned. Two files will be created in the path specified in out_file:
#		A csv file, with 5 columns:
#			nterm	cterm	AAseq		ucsd_format        prob   method
#		e.g.	AD	LEWMD	ADSLEWMD	A D S L E W M D     1     new_addin_1
#		AAseq is the full sequence of amino-acids of the peptide, prob is the predicted probability of being hit.
writePep <- function(newPep, S.Pos, pred_prob, classlist, out_file, method_name) 
{
	#Generate the amino-acid on each position 
	nPep <- dim(newPep)[1]
	nF <- dim(newPep)[2]
	pepTable <- c()
	for(i in 1:nPep) {
		nTerm <- c()
		for(j in 1:S.Pos) {
			if( newPep[i,j] != -1 ) {
				AApool <- colnames(classlist)[which(classlist == newPep[i,j])]
				AA <- AApool[sample.int(length(AApool), size = 1)]
				nTerm <- paste(nTerm, AA, sep = "")
			}
		}
		cTerm <- c()
		for(j in (S.Pos+1):nF) {
			if( newPep[i,j] != -1 ) {
				AApool <- colnames(classlist)[which(classlist == newPep[i,j])]
				AA <- AApool[sample.int(length(AApool), size = 1)]
				cTerm <- paste(cTerm, AA, sep = "")
			}
		}
		AAseq <- paste(nTerm, 'S', cTerm, sep = "")
        ucsd_seq <- paste(unlist(strsplit(AAseq, '')), collapse=' ')
		pepTable <- rbind(pepTable, c(nTerm, cTerm, AAseq, ucsd_seq))
	}
	pepTable <- cbind(pepTable, pred_prob, rep(method_name,nPep))
	colnames(pepTable) <- c('nterm', 'cterm', 'AAseq', 'ucsd format', 'prob', 'method')
	#Sort peptides according to their probabilities of being hits
	pepTable <- pepTable[order(pred_prob, decreasing = TRUE),]
	rownames(pepTable) <- c(1:nPep)
	pepTable <- as.data.frame(pepTable, stringsAsFactors = FALSE)	
	write.csv(pepTable, out_file)
}
