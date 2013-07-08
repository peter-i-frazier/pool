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