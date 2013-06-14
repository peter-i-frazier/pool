NB.predict <- function(newdata, theta, prior.positive = 10**(-4))
{
#================================================================================
#Function: NB.predict
#
#--------------------------------------------------------------------------------
#Description:
#    Predict the outcome value of each peptide in newdata.
#
#--------------------------------------------------------------------------------
#Input arguments:
#    newdata
#            A data frame whose rows correspond to peptides and columns to featu-
#            -res. The last three columns are the outcome values(whether each pe-
#            -ptide works with enzyme 1, 2a and 2b).
#    theta
#            Likelihood parameters. A list returned by the funciton getTheta.
#    prior.positive
#            Prior distribution Pr(Y=1) where Y is the outcome value. The default 
#            10^-4. Considering all peptides with length 11
#--------------------------------------------------------------------------------
#Return objects:
#    predict
#            A vector of posterior outcomes of each peptide.  
#--------------------------------------------------------------------------------
    theta_0 <- theta$theta_0
    theta_1 <- theta$theta_1
    nData <- dim(newdata)[1]
    #K: number of features
    K <- dim(newdata)[2] - 3
    predict <- rep(0, nData)
    for(n in 1:nData) {
        feature <- newdata[n,1:K]
        likelihood_1 <- 1
        likelihood_0 <- 1
        for(i in 1:K) {
            if(!is.na(feature[i])){
                likelihood_1 <- theta_1[i,feature[i]]*likelihood_1
                likelihood_0 <- theta_0[i,feature[i]]*likelihood_0
            }
        }
        predict[n] = likelihood_1*prior.positive/(likelihood_0*(1-prior.positive)
            + likelihood_1*prior.positive)
    }
    return(predict)
}
