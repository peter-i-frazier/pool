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
