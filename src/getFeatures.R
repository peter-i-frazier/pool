getFeatures <- function(data.org, AAclass, nL, nR)
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
#             res.
#    nR
#             Number of amino acids at the right of the serine that will be feat-
#             ures.
#
#--------------------------------------------------------------------------------
#Return objects:
#    feature
#             A data frame whose rows correspond to peptides and columns to featu
#             -res. The last three columns are the outcome values(whether each pe
#             -ptide works with enzyme 1, 2a and 2b).
#
#--------------------------------------------------------------------------------
    #nVal: number of values each feature can take
    nVal <- max(AAclass)
    #nOUTCOME: number of outcome values
    nOUTCOME <- dim(data.org)[2]-2
    feature <- matrix(NA, nrow=dim(data.org)[1], ncol=nL+nR+nOUTCOME)
    for (r in 1:dim(data.org)[1]) {
        sequence <- unlist(strsplit(data.org[r, 'nterm'],''))
        l.seq  <- length(sequence)
        l <- min(nL,l.seq)
        for (i in 1:l) {
		feature[r,nL+1-i] <- AAclass[1,sequence[l.seq+1-i]]
	  }
        sequence <- unlist(strsplit(data.org[r, 'cterm'],''))
        l.seq  <- length(sequence)
        c <- min(nR,l.seq)
        for (i in 1:c) {
		feature[r,nL+i] <- AAclass[1,sequence[i]]
	  }
        if( nOUTCOME != 0) {
            for (i in 1:nOUTCOME) {
                feature[r,nL+nR+i] <- data.org[r,i]
            }
        }
    }    
    #outcome.names : name of outcome values
    outcome.names <- c()
    if(nOUTCOME != 0) {
        outcome.names <- colnames(data.org)[1:nOUTCOME]
    }
    feature <- as.data.frame(feature)
    if(nOUTCOME != 0){
        colnames(feature) <- c(paste('L',nL:1,sep=""),
            paste('R',1:nR,sep=""), outcome.names)
    }
    else {
        colnames(feature) <- c(paste('L',nL:1,sep=""), 
            paste('R',1:nR,sep=""))
    }
    return( feature )
}
      
