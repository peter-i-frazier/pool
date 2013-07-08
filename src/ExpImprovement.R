ExpImprovement <- function(train, newdata, nrep, nVal, best.length = 11,
 prior.positive = 10**-4)
{
#================================================================================
#Function: ProbImprovement
#
#--------------------------------------------------------------------------------
#Description:
#    Given a set of peptides(newdata), computes the expected reduction in length 
#    of the shortest peptide known compared with current shortest length
#    (best.length).
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
#    exp.improve
#            The desired probability.                 
#
#--------------------------------------------------------------------------------
    nData <- dim(newdata)[1]
    #prob.hit: A binary vector. prob.hit[i] is the probability the ith peptide in 
	#          newdata is a hit(according to simulated posterior theta).
    prob.hit <- rep(0, nData)
	#is.hit: A binary vector. is.hit[i] = 1 if the ith peptide in newdata is a 
	#        hit. This vector is simulated according to  prob.hit.
	is.hit <- rep(0, nData)
	#reduction: A vector. reduction[i] is the reduction in the length compared
    #        	with best.length of the ith peptide in newdata if that peptide 
	#           is a hit
    #                     reduction[i] = 0 if is.hit[i] = 0
	#           reduction[i] = max(0,best.length-length(newdata[i])) 
	#           if is.hit[i] = 1
	
	#Calculate the length of all peptides in newdata
	peptide.length <- rowSums(!is.na(newdata))
	#Initialize reduction
	reduction <- pmax(0,best.length - peptide.length)
    #max.reduction: A vector. max.reduction[i] = is the maximum reduction in the 
	#               length among all peptides in newdata in the ith simulation.
    max.reduction <- rep(0, nrep)
	alpha <- Dirichlet_Parameters(train, nVal)
    for(i in 1:nrep) {        
        theta <- getTheta_MC(alpha = alpha, nVal = nVal)
		prob.hit <- NB_predict(newdata, theta, prior.positive = prior.positive)
        #Simulate is.hit according to prob.hit
		is.hit <- rbern(nData, prob.hit)
		reduction <- reduction*is.hit
		max.reduction[i] <- max(is.hit)
    }
    exp.improve <- mean(max.reduction)
    return(exp.improve)   
}
    