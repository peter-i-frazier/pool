getFeatures <- function(data.org, class, nL, nR)
{
#================================================================================
#Function: getFeatures
#
#--------------------------------------------------------------------------------
#Description:    
#    Get the feature vectors of all peptides. 
#
#--------------------------------------------------------------------------------
#Input arguments:
#    data.org
#             A matrix whose rows correspond to peptides. The matrix has five co-
#             -lumns. The first 3 shows whether each peptide works with enzyme 1, 
#             2a and 2b. The 3rd column contains amino-acids at the left of the 
#             serine and the 4th contains those at right.
#
#    class
#             A matrix specifying the mapping between each amino-acids and the 
#             class it belongs to
#    nL
#             Number of amino acids at the left of the serine that will be featu-
#             -res.
#    nR
#             Number of amino acids at the right of the serine that will be feat-
#             -ures.
#
#---------------------------------------------------------------------------------
#Return objects:
#    feature
#             A list containing an integer nVal and a data frame feature.data.
#        nVal
#             Number of values each feature can take.
#        feature.data
#             A data frame whose rows correspond to peptides and columns to featu-
#             -res. The last three columns are the outcome values(whether each pe-
#             -ptide works with enzyme 1, 2a and 2b).
#
#---------------------------------------------------------------------------------
    nVal <- max(class)
    feature.data <- matrix(NA, nrow=dim(data.org)[1], ncol=nL+nR+3)
    for (r in 1:dim(data.org)[1]) {
        sequence <- unlist(strsplit(data.org[r, 'nterm'],''))
        l.seq  <- length(sequence)
        l <- min(nL,l.seq)
        for (i in 1:l) {
		feature.data[r,nL+1-i] <- AAclass[1,sequence[l.seq+1-i]]
	  }
        sequence <- unlist(strsplit(data.org[r, 'cterm'],''))
        l.seq  <- length(sequence)
        c <- min(nR,l.seq)
        for (i in 1:c) {
		feature.data[r,nL+i] <- AAclass[1,sequence[i]]
	  }
        feature.data[r,nL+nR+1] <- data.org[r,'sfp']
        feature.data[r,nL+nR+2] <- data.org[r,'PaAcpH']
        feature.data[r,nL+nR+3] <- data.org[r,'Enzyme_2b']
    }
    feature.data <- as.data.frame(feature.data)
    colnames(feature.data) <- c(paste('L',nL:1,sep=""),paste('R',1:nR,sep=""),
           'sfp','PaAcpH','Enzyme_2b')
    feature <- list(nVal = nVal, feature.data = feature.data)
    return( feature )
}
      
