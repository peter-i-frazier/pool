Dirichlet_Parameter <- function(train, nVal) 
{
#===============================================================================
#Function: Dirichlet_Parameter
#-------------------------------------------------------------------------------
#Description:
#    Computes the posterior parameter vector
#                       alpha = (alpha_1,alpha_2,...alpha_K) of 
#    Dirichlet distribution Dir(alpha,K) for each feature given training data. 
#    The prior parameter vector of a feature d positions away from the serine is
#    set as
#                       (d**0.5,...,d**0.5)
# 
#-------------------------------------------------------------------------------
#Input arguments:
#    train
#             A data matrix. Each row corresponds to a peptide and each column 
#             correpsonds to a feature.
#    nVal
#             Number of values each feature can take.
#-------------------------------------------------------------------------------
#Return objects:
#    alpha
#             A data matrix. The rows correspond to possible values a feature 
#             can take and columns to features. Each entry (x,y) 
#             corresponds to the alpha parameter that feature y takes value 
#             x.
#
#-------------------------------------------------------------------------------
      # K: number of features
      K <- dim(train)[2] 
           
	R <- dim(train)[1]
      alpha <- c()
	for (col in 1:K) {
		count <- rep(0,nVal)
		for (r in 1:R) {
                  if (!is.na(train[r,col])) {
                      count[train[r,col]] <- count[train[r,col]] + 1
                  }
		}
            alpha <- cbind(alpha,count + rep(abs(col-(K+1)/2)**0.5,nVal))
      }	
      return (alpha)
}