get_P <- function(prob_vec) {
    num_rec <- length(prob_vec)
    PP <- rep(0, num_rec)
    prod <- 1
    for (n in 1:num_rec) {
        prod <- prod * (1 - prob_vec[n])
        PP[n] <- 1 - prod
    }
    return (PP)
}
