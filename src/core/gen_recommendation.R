# Author: Jialei Wang
# Created: 05.29.2015
# Functions for generating recommended peptides using various methods

GenOnePeptideMAPOld <- function(length.left, length.right, alpha.1, alpha.0, num.mc.samples) {
  # Construct feature vector for one peptide sequence using Maximize A Posteriori,
  # the way to construct it is that we use estimated posterior distribution over 
  # thetas, and for each position, we pick the AA with highest theta.1 / theta.0
  #
  # Args:
  #   length.left: int, length of the sequence to be constructed to the left of Serine.
  #   length.right: int, length of the sequence to be constructed to the right of Serine.
  #   alpha.1: returned by BayesianNaiveBayes, parameter of the posterior Dirichlet distribution.
  #   alpha.0: returned by BayesianNaiveBayes, parameter of the posterior Dirichlet distribution.
  #   num.mc.samples: int, number of samples of thetas sampled from posterior distribution.
  #
  # Returns:
  #   vector, feature vector of the constructed peptide
  feature <- rep(-1, MAXL + MAXR)
  thetas.1 <- SampleThetas(alpha.1, num.mc.samples)
  thetas.0 <- SampleThetas(alpha.0, num.mc.samples)
  for (j in c(1:length.left, MAXL + 1:length.right)) {
    ratio <- thetas.1[[j]] / thetas.0[[j]]
    index.max <- sapply(1:num.mc.samples, function (i) which.max(ratio[i, ]))
    # find mode of index.max and set feature[j] to it
    count.of.index.max <- sapply(1:ncol(ratio), function (i) sum(index.max == i))
    feature[j] <- which.max(count.of.index.max)
  }
  return (feature)
}

GenOnePeptideMAPNew <- function(length.left, length.right, alpha.1.label, alpha.0.label,
                                alpha.1.not.label, alpha.0.not.label, alpha.1.unlabel, 
                                alpha.0.unlabel, num.mc.samples) {
  # The only difference from GenOnePeptideMAPOld is that, we train three classifiers,
  # each predicting one step in the experiment: label, not label and unlabel.
  #
  # Args:
  #   length.left: int, length of the sequence to be constructed to the left of Serine.
  #   length.right: int, length of the sequence to be constructed to the right of Serine.
  #   alpha.1.label: returned by BayesianNaiveBayes, for predicting "label" step.
  #   alpha.0.label: returned by BayesianNaiveBayes, for predicting "label" step.
  #   alpha.1.not.label: returned by BayesianNaiveBayes, for predicting "not label" step.
  #   alpha.0.not.label: returned by BayesianNaiveBayes, for predicting "not label" step.
  #   alpha.1.unlabel: returned by BayesianNaiveBayes, for predicting "unlabel" step.
  #   alpha.0.unlabel: returned by BayesianNaiveBayes, for predicting "unlabel" step.
  #   num.mc.samples: int, number of samples of thetas sampled from posterior distribution.
  #
  # Returns:
  #   vector, feature vector of the constructed peptide
  feature <- rep(-1, MAXL + MAXR)
  thetas.1.label <- SampleThetas(alpha.1.label, num.mc.samples)
  thetas.0.label <- SampleThetas(alpha.0.label, num.mc.samples)
  thetas.1.not.label <- SampleThetas(alpha.1.not.label, num.mc.samples)
  thetas.0.not.label <- SampleThetas(alpha.0.not.label, num.mc.samples)
  thetas.1.unlabel <- SampleThetas(alpha.1.unlabel, num.mc.samples)
  thetas.0.unlabel <- SampleThetas(alpha.0.unlabel, num.mc.samples)
  for (j in c(1:length.left, MAXL + 1:length.right)) {
    ratio <- thetas.1.label[[j]] / thetas.0.label[[j]] * thetas.0.not.label[[j]] / 
             thetas.1.not.label[[j]] * thetas.1.unlabel[[j]] / thetas.0.unlabel[[j]]
    index.max <- sapply(1:num.mc.samples, function (i) which.max(ratio[i, ]))
    # find mode of index.max and set feature[j] to it
    count.of.index.max <- sapply(1:ncol(ratio), function (i) sum(index.max == i))
    feature[j] <- which.max(count.of.index.max)
  }
  return (feature)
}

GenRecomSetOld <- function(X, Y, alpha.1, alpha.0, p1, num.recom, num.mc.samples, 
                           minL, maxL, minR, maxR)  {
  # Generates a recommendation set
  # 
  # Args:
  #   X: Feature matrix of peptides.
  #   Y: Label vectors of peptides.
  #   alpha.1: List of vectors, where each vector is parameter for a Dirichlet prior
  #            corresponding to one feature for x with y = 1.
  #   alpha.0: List of vectors, where each vector is parameter for a Dirichlet prior
  #            corresponding to one feature for x with y = 0.
  #   p1: P(y = 1), which is a constant
  #   num.recom: Number of recommendations to generate.
  #   num.mc.samples: Number samples to draw.
  #   minL: Min length allowed for the left half of a peptide.
  #   maxL: Max length allowed for the left half of a peptide.
  #   minR: Min length allowed for the right half of a peptide.
  #   maxR: Max length allowed for the right half of a peptide.
  #
  # Returns:
  #   A matrix of recommendation set.
  if (num.recom < 1 || num.recom %% 1 > 0)
    stop("num.recom must be a positive integer")
  num.unique.recom <- 0
  recom.set <- c()
  while (num.unique.recom < num.recom) { 
    print(sprintf("iter %d", num.unique.recom))

    alphas <- BayesianNaiveBayes(X, Y, alpha.1, alpha.0, p1)
    length.left <- ceiling(runif(1, min = minL - 1, max = maxL))
    length.right <- ceiling(runif(1, min = minR - 1, max = maxR))
    new.rec.peptide <- GenOnePeptideMAPOld(length.left, length.right,
                                           alphas$post.alpha.1,
                                           alphas$post.alpha.0, num.mc.samples)
    if (nrow(unique(rbind(recom.set, new.rec.peptide))) ==
        num.unique.recom)  # already contained this peptide
      next
    num.unique.recom <- num.unique.recom + 1
    X <- rbind(X, new.rec.peptide)
    Y <- c(Y, 0)
    recom.set <- rbind(recom.set, new.rec.peptide)
  }
  return (recom.set)
}

GenOneRandomPeptide <- function(minL, maxL, minR, maxR, max.left,
                                max.right, class.vec) {
  # Generates a random peptide
  # 
  # Args:
  #   minL: Min length allowed for the left half of a peptide.
  #   maxL: Max length allowed for the left half of a peptide.
  #   minR: Min length allowed for the right half of a peptide.
  #   maxR: Max length allowed for the right half of a peptide.
  #   max.left: Hard limit on the length of the left half of a peptide
  #   max.right: Hard limit on the length of the right half of a peptide
  #              i.e. the total length of a peptide is (max.left + max.right)
  #   class.vec: Vector of class numbers of amino acids
  #
  # Returns:
  #   Feature vector of a peptide
  length.left <- sample(minL:maxL, 1)
  length.right <- sample(minR:maxR, 1)
  feature <- rep(-1, max.left + max.right)
  feature[1:length.left] <- sample(class.vec, length.left, replace = TRUE)
  feature[max.left + 1:length.right] <- sample(class.vec, length.right,
                                               replace = TRUE)
  return (feature)
}

