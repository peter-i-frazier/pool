# Author: Jialei Wang
# Created: 05.28.2015
# This file contains all functions related to modeling reversible labeling project,
# including retrieving data, and setting up training model.

GetDataReducedAA <- function(query) {
  # Retrieve data from db, and parse information, including feature extraction.
  # 
  # Args:
  #   query: SQL query of retrieving data
  #
  # Returns:
  #   list containing feature table and data retrieved from query
  data <- QuerySql("sfp_AcpS", query)
  data <- data[data[, 'nterm'] != "", ] # filter out rows that do not have nterm & cterm
  for (r in 1:nrow(data)) {
    if (paste0(data[r, 'nterm'], "S", data[r, 'cterm']) != data[r, 'sequence'])
      print ("error: inconsistent sequence! reversible_labeling_model.R L18")
  }
  feature <- GetFeatureFromSeqsReducedAA(data[, 'nterm'], data[, 'cterm'])
  return (list(feature=feature, data=data))
}

GetFeatureFromSeqsReducedAA <- function(nterms, cterms) {
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
  #
  # Returns:
  #   Matrix, with each row corresponding to feature vector of one sequence.
  lookup.table <- QuerySql("sfp_AcpS", "SELECT * FROM reduced_AA")
  if (length(nterms) != length(cterms))
    print ("error: length of nterms and cterms are not equal!")
  feature.table <- matrix(-1, nrow=length(nterms), ncol=38)
  nterms.list <- strsplit(nterms, split='')
  cterms.list <- strsplit(cterms, split='')
  for (n in 1:length(nterms)) {
    len.nterm <- length(nterms.list[[n]])
    len.cterm <- length(cterms.list[[n]])
    for (j in 1:min(MAXL, len.nterm))
      feature.table[n, j] <- lookup.table[lookup.table[, 'AA'] == nterms.list[[n]][len.nterm - j + 1], 'class']
    for (j in 1:min(MAXR, len.cterm))
      feature.table[n, j+MAXL] <- lookup.table[lookup.table[, 'AA'] == cterms.list[[n]][j], 'class']
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
