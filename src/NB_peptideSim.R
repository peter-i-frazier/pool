NB_peptideSim <- function(train = NA, theta = NA, outcome_name = NA, nVal, nPep, 
isHit = 1, minL = 1, maxL = 3, minR = 3, maxR = 9 )
{
#================================================================================
#Function: NB_predict
#
#--------------------------------------------------------------------------------
#Description:
#         Simulate peptides according to the posterior distribution of theta. 
#    Theta can be specified or trained using training data. 
#    The length of simulated peptides:
#        Currently, the length to the left & the right of the serene are drawn 
#    from independent uniform distributions, where the left is uniform from 
#    length minL(default value is 1) to maxL(default value is 3), and the right is 
#    uniform from length minR(default value is 3) to maxR(default value is 9).
#        For further improvement, the distribution of the length of peptides 
#    could be a parameter of this function.
#
#--------------------------------------------------------------------------------
#Input arguments:
#    train
#             A data matrix. Each row corresponds to a peptide and each column 
#             corresponds to a feature with the last column being the outcome
#             value.
#    theta
#             A list of two likelihood matrices: theta_0 and theta_1. 
#             The rows of a likelihood matrix correspond to possible values a 
#             feature can take and columns to features. Each entry (x,y) 
#             corresponds to a likelihood probability that feature y takes value 
#             x. 
#    The user may either specify train or theta, or neither. If both are 
#    specified, TRAIN is ignored. If neither is specified, peptides will be
#    outcome
#             Name of the outcome value.  
#    nVal
#             Number of values each feature can take.
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
peptides <- matrix("",nrow = nPep, ncol = 3)
if(missing(theta) && !missing(train)){
	#Get alpha parameter from training data
    alpha <- Dirichlet_Parameter(train,nVal)
}
for(j in 1:nPep) {
    #Draw the left & right lengths at random from the uniform distribution\
    sL <- minL-1+ceiling((maxL-minL+1)*runif(1))
    sR <- minR-1+ceiling((maxR-minR+1)*runif(1))
	if(!(missing(theta) && missing(train))) {
	    if(missing(theta)){
			#Simulate theta according to alpha parameter
			theta <- getTheta_MC(alpha = alpha, nVal = nVal)	
	    }
        if(isHit) {
	        param <- theta$theta_1 }
	    else {
			param <- theta$theta_0 }		
	}	
    #Iterate through positions
    for(i in sL:1) {
	    if(missing(theta) && missing(train)) {
			s_class <- ceiling(8*runif(1)) 
		}
		else {
			s_class <-  rdist(c(1:nVal), param[,paste('L',i,sep="")])
		}
        s_class_list <- colnames(class)[class == s_class]
        s_amino_acid <- s_class_list[ceiling(length(s_class_list)*runif(1))]
        peptides[j,2] <- paste(peptides[j,2], s_amino_acid, sep="")
		}
    for(i in 1:sR) {
		if(missing(theta) && missing(train)) {
			s_class <- ceiling(8*runif(1)) 
		}
		else {
			s_class <-  rdist(c(1:nVal), param[,paste('R',i,sep="")])
		}
        s_class_list <- colnames(class)[class == s_class]
        s_amino_acid <- s_class_list[ceiling(length(s_class_list)*runif(1))]
        peptides[j,3] <- paste(peptides[j,3], s_amino_acid, sep="")
    }    
}
peptides[,1] <- rep(isHit, nPep)
colnames(peptides) <- c(outcome_name,'nterm','cterm')
return(peptides)
}
