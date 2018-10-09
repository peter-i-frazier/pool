# Author: Jialei Wang
# Created: 05.28.2015
# This file contains all functions related to modeling reversible labeling project,
# including retrieving data, and setting up training model.

GetDataReducedAA <- function(table, lookup) {
  # parse table, including feature extraction.
  # 
  # Args:
  #   table: data table
  #
  # Returns:
  #   list containing feature table and data retrieved from query
  data <- data[data[, 'nterm'] != "", ] # filter out rows that do not have nterm & cterm
  feature <- GetFeatureFromSeqsReducedAA(as.character(data[, 'nterm']), as.character(data[, 'cterm']), lookup)
  return (list(feature=feature, data=data))
}

GetFeatureFromSeqsReducedAA <- function(nterms, cterms, lookup_table) {
  # Extract feature vectors from sequences, where each sequence contains one nterm
  # and cterm. The vector is of length (MAXL + MAXR), with (1 : MAXL)th entries 
  # extracted from nterm and (MAXL+1 : MAXL+MAXR)th entries extracted from cterm.
  # Note: for nterm, we extract feature from right to left, because in designing
  # prior, we consider dependency upon distance of the AA to the Serine!
  #
  # Args:
  #   nterms: vector of strings, where each element corresponds to N-terminus of a
  #           sequence.
  #   cterms: vector of strings, corresponds to C-terminus
  #   lookup_table: lookup between AA and class
  #
  # Returns:
  #   Matrix, with each row corresponding to feature vector of one sequence.
  if (length(nterms) != length(cterms))
    print ("error: length of nterms and cterms are not equal!")
  feature.table <- matrix(-1, nrow=length(nterms), ncol=38)
  nterms.list <- strsplit(nterms, split='')
  cterms.list <- strsplit(cterms, split='')
  for (n in 1:length(nterms)) {
    len.nterm <- length(nterms.list[[n]])
    len.cterm <- length(cterms.list[[n]])
    for (j in 1:min(MAXL, len.nterm))
      feature.table[n, j] <- lookup_table[lookup_table[, 'AA'] == nterms.list[[n]][len.nterm - j + 1], 'class']
    for (j in 1:min(MAXR, len.cterm))
      feature.table[n, j+MAXL] <- lookup_table[lookup_table[, 'AA'] == cterms.list[[n]][j], 'class']
  }
  return (feature.table)
}

SetPriorReducedAA <- function(coefficient, num.class) {
  # Calculate parameters for a prior distribution, and it is designed such that
  # alpha is getting bigger when the position is away from Serine.
  #
  # Args:
  #   coefficient: double, set the scale of the alphas.
  #   num.class: int, number of classes for reducedAA
  #
  # Returns a list:
  #   jth element: vector, is alpha for jth position.
  alpha <- list()
  for (j in 1:MAXL)
    alpha[[j]] <- coefficient * rep(sqrt(j), num.class)
  for (j in (MAXL+1):(MAXL+MAXR))
    alpha[[j]] <- coefficient * rep(sqrt(j - MAXL), num.class)
  return (alpha)
}

GenRandomSeqFromFeature <- function(feature, lookup.table) {
  # Given feature of a peptide, generate the amino acid sequences, where Serine
  # should not be in the sequence except for the anchor Serine.
  # return a list with nterm and cterm of the generated peptide
  nterm <- c()
  for (i in 1:MAXL) {
    if (feature[i] == 7) {
      nterm <- c("T", nterm)
    } else if (feature[i] > 0) {
      nterm <- c(sample(lookup.table[lookup.table[, 'class'] == feature[i], 'AA'], 1), nterm)
    }
  }
  cterm <- c()
  for (i in 1:MAXR) {
    if (feature[i + MAXL] == 7) {
      cterm <- c(cterm, "T")
    } else if (feature[i + MAXL] > 0) {
      cterm <- c(cterm, sample(lookup.table[lookup.table[, 'class'] == feature[i + MAXL], 'AA'], 1))
    }
  }
  return (list(nterm = paste0(nterm, collapse=''), cterm = paste0(cterm, collapse='')))
}
