getTheta <- function(train) {
	K <- dim(train)[2]
	R <- dim(train)[1]
	theta <- c()

	for (col in 1:K) {
		count <- rep(0,6)
		for (r in 1:R) {
			count[train[r,col]] <- count[train[r,col]] + 1
			alpha <- count + rep(abs(col-19.5)**0.5,6)
		}
			theta <- cbind(theta, alpha/sum(alpha))
	}
	return (theta)
}