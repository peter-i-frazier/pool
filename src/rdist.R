rdist <- function(value, p) 
{
#===============================================================================
#Function: rdist
#
#-------------------------------------------------------------------------------
#Description:
#   Generate a random number that has a discrete distribution.
# 
#-------------------------------------------------------------------------------
##Input arguments:
#   value
#             A vector contain the possiblly taken values of the random number.
#   p
#             A vector contain the probability the random number taking each 
#             value. p must have a 1-1 correspondence with value.
#-------------------------------------------------------------------------------
#Return objects:
#   rdn
#             Generated random number.             
#
#-------------------------------------------------------------------------------
    l <- length(p)
    sum_p <- rep(0, l)
    for(i in 1:l) {
        sum_p[i] <- sum(p[1:i])
    }
    u_rdn <- runif(1)
    for(i in 1:l) {
        if(u_rdn < sum_p[i]){
            rdn <- value[i]
            break
        }
    }
    return(rdn)
}
