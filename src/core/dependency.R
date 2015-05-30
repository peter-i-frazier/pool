# Author: Jialei Wang
# Created: 05.26.2015
# Resolve dependency

ResolveDependency <- function(root) {
  source(paste0(root, "core/constants.R"))
  source(paste0(root, "core/util.R"))
  source(paste0(root, "core/bayesian_naive_bayes.R"))
  source(paste0(root, "core/reversible_labeling_model.R"))
  source(paste0(root, "core/gen_recommendation.R"))
}

