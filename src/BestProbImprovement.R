ProbImprovement <- function(train, newdata, nrep)
{
#================================================================================
#Function: ProbImprovement
#
#--------------------------------------------------------------------------------
#Description:
#    Computes the probability that at least one peptide in newdata is a hit. 
#
#--------------------------------------------------------------------------------
#Imput arguments:
#    train:
#            A list returned by the function getFeatures.
#    newdata:
#            A data frame whose rows correspond to peptides and columns to featu-
#            -res. The last three columns are the outcome values(whether each pe-
#            -ptide works with enzyme 1, 2a and 2b).
#    nrep:
#            Number of simulation replications.
#--------------------------------------------------------------------------------
#Return objects:
#    prob.improve
#            The desired probability.                 
#
#--------------------------------------------------------------------------------
    trainData <- train$feature.data
    nVal <- train$nVal
    #train_X: A data matrix with a feature vector at each row and the outcome va-
    #         -lue X at the last column. X = 1, 2a, 2b 
    train_1 <- trainData[,1:dim(trainData)[2]-2]
    train_2a <- trainData[,c(1:dim(trainData)[2]-3,dim(trainData)[2]-1)]
    train_2b <- trainData[,c(1:dim(trainData)[2]-3,dim(trainData)[2])]

    #prob.hit: A vector. prob.hit[i] is the probability the ith peptide in newda-
    #          -is a hit.
    prob.hit <- rep(0, dim(newdata)[1])
    #prob.find: The probability that among all peptides in newdata at least one 
    #           is a hit.
    prob.find <- rep(0, nrep)
    for(i in 1:nrep) {        
        theta_1 <- getTheta.MC(train_1, nVal)
        predict_1 <- NB.predict(newdata, theta_1)        
        theta_2a <- getTheta.MC(train_2a, nVal)
        predict_2a <- NB.predict(newdata, theta_2a)        
        theta_2b <- getTheta.MC(train_2b, nVal)
        predict_2b <- NB.predict(newdata, theta_2b)
        for(x in 1:dim(newdata)[1]){
            prob.hit[x] <- predict_1[x]*(1-(1-predict_2a[x])*(1-predict_2b[x]))
        }
        prob.find[i] <- 1-prod(1-prob.hit)
    }
    prob.improve <- sum(prob.find)/nrep
    return(prob.improve)   
}
    
    