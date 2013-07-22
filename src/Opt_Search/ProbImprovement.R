ProbImprovement <- function(train, newdata, nrep, nVal, best.length = 11,
 prior.positive = 10**-4)
{
#================================================================================
#Function: ProbImprovement
#
#--------------------------------------------------------------------------------
#Description:
#    Computes the probability that at least one peptide with length shorter than 
#    previous best in newdata is a hit.
#
#--------------------------------------------------------------------------------
#Input arguments:
#    train
#            A data matrix. Each row corresponds to a peptide and each column 
#            corresponds to a feature with the last column being the outcome
#            value.
#    newdata:
#            A data matrix whose rows correspond to peptides and columns to 
#            features. There shall be no outcome values.
#    nrep:
#            Number of simulation replications.
#    nVal:
#            Number of values each feature can take.
#    best.length
#            Current shortest length of peptides working with both enzymes.
#    prior.positive
#            Prior distribution Pr(Y=1) where Y is the outcome value. The default 
#            10^-4. 
#
#--------------------------------------------------------------------------------
#Return objects:
#    prob
#        A list which contains:
#    prob.improve
#            The desired expected reduction. 
#    prod.CI                
#            The 95% confidence interval of the desired probability.
#--------------------------------------------------------------------------------
    nData <- dim(newdata)[1]
    #prob.hit: A binary vector. prob.hit[i] is the probability the ith peptide in 
	#          newdata is a hit(according to simulated posterior theta).
    prob.hit <- rep(0, nData)
	#is.hit: A binary vector. is.hit[i] = 1 if the ith peptide in newdata is a 
	#        hit, this vector is simulated according to  prob.hit
	is.hit <- rep(0, nData)
	#isShorter: A vector. isShorter[i] is the indicator of the ith peptide in 
	#           newdata being shorter than the best previous peptide.
	is.Shorter <- as.numeric(newdata < best.length)
    #is.found: A binary vector. In the ith simulation, if at least one of the 
	#          peptides in newdata is a hit and is shorter than best.length, 
	#      	   is.found[i] = 1, otherwise 0.
    is.found <- rep(0, nrep)
	alpha <- Dirichlet_Parameter(train, nVal)
    for(i in 1:nrep) {        
        theta <- getTheta_MC(alpha = alpha, nVal = nVal)
		prob.hit <- NB_predict(newdata, theta, prior.positive = prior.positive)
        #Simulate is.hit according to prob.hit
		is.hit <- rbern(nData, prob.hit)
		is.found[i] <- as.numeric(any(is.hit*is.shorter > 0))
    }
    prob.improve <- mean(is.found)
	prob.CI <- c(prob.improve-sd(is.found)*qnorm(0.975),prob.improve+sd(is.found)*qnorm(0.975))
	prob <- list("prob.improve" = prob.improve, "prob.CI" = prob.CI)
    return(prob)   
}
    