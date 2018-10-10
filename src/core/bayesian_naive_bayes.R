# Author: Jialei Wang
# Created: 05.28.2015
# This file provides utility functions for Bayesian Naive Bayes classifier. 
# The model is stated in the following:
#
# For any x, where x is a vector, and its elements is drawn independently
# from some categorical distributions which are independent across position,
# parameterized by \theta_j, and j is index for position. If we assume
# \theta_j has prior distribution which is Dirichlet, then posterior distribution
# of \theta_j given x is also Dirichlet. Now suppose each x has a label y \in {0, 1}, 
# and let \theta_j^y determine the categorical distributions that generate x having
# label y, then using Bayes' formula, we can calculate
# P(y = 1 | x)
# = P(x | y = 1) * P(y = 1) / ( P(x | y = 1) * P(y = 1) + P(x | y = 0) * P(y = 0) )
# = \Prod_j \theta_j^1 * P(y = 1) / ( \Prod_j \theta_j^1 * P(y = 1) + \Prod_j \theta_j^0 * P(y = 0) )
# Notationwise, let \theta_j^y have prior distribution Dirichlet(\alpha_j^y)

require(MCMCpack)

BayesianNaiveBayes <- function(X, Y, alpha.1, alpha.0, p1) {
  # Construct a Bayesian Naive Bayes model.
  #
  # Args:
  #   X: matrix with each row corresponding to feature vector of one x.
  #   Y: vector of labels {0, 1}
  #   alpha.1: list of vectors, where each vector is parameter for a Dirichlet prior
  #            corresponding to one feature for x with y = 1.
  #   alpha.0: list of vectors, where each vector is parameter for a Dirichlet prior
  #            corresponding to one feature for x with y = 0.
  #   p1: P(y = 1), which is a constant
  # Returns:
  #   list that contains post.alpha.1 (list), post.alpha.0 (list), p1 (double)
  X.1 <- X[Y == 1,]
  if (is.null(nrow(X.1))) {
    X.1 <- matrix(X.1, nrow=1, ncol=length(X.1))
  }
  X.0 <- X[Y == 0,]
  if (is.null(nrow(X.0))) {
    X.0 <- matrix(X.0, nrow=1, ncol=length(X.0))
  }
  if (length(alpha.1) != length(alpha.0))
    stop ("error: length(alpha) are not consistent!")
  for (j in 1:length(alpha.1)) {
    if (length(alpha.1[[j]]) != length(alpha.0[[j]]))
      stop ("error: length(alpha[[j]]) are not consistent!")
    for (i in 1:length(alpha.1[[j]])) {
      alpha.1[[j]][i] <- alpha.1[[j]][i] + sum(X.1[, j] == i)
      alpha.0[[j]][i] <- alpha.0[[j]][i] + sum(X.0[, j] == i)
    }
  }
  return (list(post.alpha.1=alpha.1, post.alpha.0=alpha.0, p1=p1))
}

SampleThetas <- function(alpha, num.samples) {
  # Draw samples of theta from Dirichlet distribution parameterized by alpha
  #
  # Args:
  #   alpha: list, with each component parameter that specifying a Dirichlet distribution.
  #   num.samples: number of samples to draw
  #
  # Returns a list:
  #   jth element: matrix, has num.samples rows, and each row is a sample of theta drawn
  #                from Dirichlet distribution
  thetas <- list()
  for (j in 1:length(alpha))
    thetas[[j]] <- rdirichlet(num.samples, alpha[[j]])
  return (thetas)
}

CalculateProb <- function(feature.matrix, thetas.1, thetas.0, p1) {
  # Based on the sample of thetas drawn from their posterior, estimate
  # distribution of P(y = 1 | x) for any x \in {x1, x2, ..., xN}
  # Note: this function runs slow when N is big, i.e., N > 1e6, and
  # should be implemented seperately for the case when N is big.
  #
  # Args:
  #   feature.matrix: matrix, each row is feature vector for one peptide
  #                   that needs to calculate P(y = 1 | x)
  #   thetas.1: list, returned by SampleThetas
  #   thetas.0: list, returned by SampleThetas
  #   p1: double, P(y=1)
  #
  # Returns a list:
  #   mean: vector, mean of estimated P(y = 1 | x) for each x
  #   sd: vector, sd of estimated P(y = 1 | x) for each x
  if (ncol(feature.matrix) != length(thetas.1) || length(thetas.1) != length(thetas.0)) {
    print(feature.matrix)
    print(thetas.1)
    stop ("error: inconsistent size!")
  }
  num.mc.samples <- nrow(thetas.1[[1]])
  num.peptides <- nrow(feature.matrix)
  p.x.1 <- matrix(1, nrow=num.mc.samples, ncol=num.peptides)
  p.x.0 <- matrix(1, nrow=num.mc.samples, ncol=num.peptides)
  for (j in 1:length(thetas.1)) {
    if (nrow(thetas.1[[j]]) != num.mc.samples || nrow(thetas.0[[j]]) != num.mc.samples)
      stop ("error: inconsistent size!")
    for (n in 1:num.peptides) {
      if (feature.matrix[n, j] != -1) {
        p.x.1[, n] <- p.x.1[, n] * thetas.1[[j]][, feature.matrix[n, j] ]
        p.x.0[, n] <- p.x.0[, n] * thetas.0[[j]][, feature.matrix[n, j] ]
      }
    }
  }
  prob <- p.x.1 * p1 / (p.x.1 * p1 + p.x.0 * (1. - p1))
  prob.mean <- sapply(1:num.peptides, function (n) mean(prob[, n]))
  prob.sd <- sapply(1:num.peptides, function (n) sd(prob[, n]))
  return (list(mean=prob.mean, sd=prob.sd))
}

Predict <- function(train_X, train_Y, test_X, alpha.1, alpha.0, p1, num.mc) {
  trained.params <- BayesianNaiveBayes(train_X, train_Y, alpha.1, alpha.0, p1)
  thetas.1 <- SampleThetas(trained.params$post.alpha.1, num.mc)
  thetas.0 <- SampleThetas(trained.params$post.alpha.0, num.mc)
  prob.list <- CalculateProb(test_X, thetas.1, thetas.0, trained.params$p1)
  return (prob.list$mean)
}


