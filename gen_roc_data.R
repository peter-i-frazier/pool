# Author: Jialei Wang
# Created: 07.22.2015
# Generate data for ROC plot

root <- "src/"
source(paste0(root, "core/dependency.R"))
ResolveDependency(root)

data <- c()
for (i in 0:4) {
  data <- rbind(data, read.csv(sprintf('data/TS%d.csv', i)))
}
AA_lookup <- read.csv('data/reduced_AA_lookup.csv')
AA_lookup[,'AA'] = as.character(AA_lookup[,'AA'])

dataset <- GetDataReducedAA(data, AA_lookup)

mc_itr = 100

CrossValidation <- function(enzyme, num_leave_out, alpha.0.prior.coef, alpha.1.prior.coef, p1, fname) {
  alpha.0.prior <- SetPriorReducedAA(alpha.0.prior.coef, NUM_CLASS)
  alpha.1.prior <- SetPriorReducedAA(alpha.1.prior.coef, NUM_CLASS)
  dataset_X <- dataset$feature[dataset$data[, enzyme] != -1, ]
  dataset_Y <- dataset$data[dataset$data[, enzyme] != -1, enzyme]
  cv_result <- leave_n_out_cv(dataset_X, dataset_Y, alpha.1.prior, alpha.0.prior, p1, mc_itr, num_leave_out)
  roc.xy <- GetRocXy(cv_result$prob, cv_result$Y)
  write.csv(roc.xy, fname, row.names=FALSE)
}

# if (TOPLOT == 1) {
#   CrossValidation('sfp', 3, 10, 0.1, 1e-5, TRUE)
#   CrossValidation('AcpS', 2, 10, 1, 1e-5, TRUE)
#   CrossValidation('PfAcpH', 3, 10, 10, 0.5, TRUE)
# } else {
#   CrossValidation('sfp', 3, ALPHA0, ALPHA1, P1, FALSE)
#   CrossValidation('AcpS', 2, ALPHA0, ALPHA1, P1, FALSE)
#   CrossValidation('PfAcpH', 3, ALPHA0, ALPHA1, P1, FALSE)
# }

print("making roc data for sfp specific")
CrossValidation('sfp_specific', 1, 10, 1, 1e-4, 'sfp_specific_roc_data.csv')
print("making roc data for acps specific")
CrossValidation('acps_specific', 1, 10, 1, 1e-5, 'acps_specific_roc_data.csv')
