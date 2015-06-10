# Author: Tom Fei
# Created: 06.06.2015
# Compute probability improvement of a set, PI(S)

ComputeProbImprovOfSet <- function(S, X, Y, alpha.1, alpha.0, p1, num.samples) {
  # Compute probability improvement of a set, S
  #
  # Args:
  #   S: A matrix of recommendation set.
  #   X: Feature matrix of peptides.
  #   Y: Label vectors of peptides.
  #   alpha.1: List of vectors, where each vector is parameter for a Dirichlet prior
  #            corresponding to one feature for x with y = 1.
  #   alpha.0: List of vectors, where each vector is parameter for a Dirichlet prior
  #            corresponding to one feature for x with y = 0.
  #   p1: P(y = 1), which is a constant
  #   num.samples: Number of samples to draw.
  # 
  # Returns:
  #   Probability improvement of the set S.
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
