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
peptides <- cbind(rep(isHit, nPep),peptides)
colnames(peptides) <- c(outcome_name,'nterm','cterm')
return(peptides)
}