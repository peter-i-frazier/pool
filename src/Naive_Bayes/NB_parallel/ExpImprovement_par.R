ExpImprovement_par <- function(train, newdata, classlist, nRep = 100, best.length = 11,
 prior.positive = 10**-4)
{
#================================================================================
#Function: ExpImprovement_par
#
#--------------------------------------------------------------------------------
#Description:
#    The parallel version of ExpImprovement. 
#        Given a set of peptides(newdata), computes the expected reduction in length 
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
#    classlist
#             A matrix specifying the mapping between each amino-acids and the 
#             class it belongs to.
#    nRep:
#            Number of simulation replications.
#    best.length
#            Current shortest length of peptides working with both enzymes.
#    prior.positive
#            Prior distribution Pr(Y=1) where Y is the outcome value. The default 
#            10^-4. 
#
#--------------------------------------------------------------------------------
#Return objects:
#	 exp_reduc:
#        A list which contains:
#    improve
#            The desired expectation of reduction in length.
#    CI                 
#            The 95% percent confidence interval of exp.improve
#
#--------------------------------------------------------------------------------
    require(Rlab)
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
	alpha <- Dirichlet_Parameter(train, classlist)
	#max_reduction: A vector. max.reduction[i] = is the maximum reduction in the 
	#               length among all peptides in newdata in the ith simulation.
    max_reduction <- foreach(icount(nRep), .init = c(), .combine = append, 
			.packages = c('Rlab'), .export = c('getTheta_MC', 'NB_predict')) %dopar% {        
        theta <- getTheta_MC(alpha = alpha, classlist = classlist)
		prob.hit <- NB_predict(newdata, theta, prior.positive = prior.positive)
        #Simulate is.hit according to prob.hit
		is.hit <- rbern(nData, prob.hit)
		max(reduction*is.hit)
    }
    improve <- mean(max_reduction)
	CI <- c(improve-sd(max_reduction)*qnorm(0.975)/sqrt(nRep), improve+sd(max_reduction)*qnorm(0.975)/sqrt(nRep))
	exp_reduc <- list("improve" = improve, "CI" = CI) 
    return(exp_reduc)   
}