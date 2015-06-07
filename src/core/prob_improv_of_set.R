# Author: Tom Fei
# Created: 06.06.2015
# Compute probability improvement of a set, PI(S)

ComputeProbImprovOfSet <- function(S, X, Y, alpha.1, alpha.0, p1, num.samples) {
  # Compute probability improvement of a set, S
  #
  # Args:
  #   S: matrix, recommendation set.
  #   The other arguments are explained in BayesianNaiveBayes
  print("Computing PI(S)")
  sum <- 0  # cumulates probabilities to get prob improv
  grow.product <- 1

  for (i in 1:nrow(S)) {
    test.x <- S[i, ]
    # ptm.0 <- proc.time()

    alphas <- BayesianNaiveBayes(X, Y, alpha.1, alpha.0, p1)
    thetas.1 <- SampleThetas(alphas$post.alpha.1, num.samples)
    thetas.0 <- SampleThetas(alphas$post.alpha.0, num.samples)

    # Note that test.x is a vector. Need to convert it to
    # "row" matrix to pass into CalculateProb.
    prob.list <- CalculateProb(t(test.x), thetas.1, thetas.0, p1)
    prob.positive <- prob.list$mean

    sum <- sum + grow.product * prob.positive
    X <- rbind(X, test.x)
    Y <- c(Y, 0)
    grow.product <- grow.product * (1 - prob.positive)

    # print(sprintf("Updating values costs %f", proc.time() - ptm.1))
  }
  return(sum)
}
