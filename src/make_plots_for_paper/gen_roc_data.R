# Author: Jialei Wang
# Created: 07.22.2015
# Generate data for ROC plot

#!/usr/local/R/icse/bin/Rscript
rm(list = ls())
root <- "/fs/home/jw865/remote_deployment/ucsd_reversible_labeling/"
source(paste0(root, "core/dependency.R"))
ResolveDependency(root)

args <- commandArgs(TRUE)
for(i in 1:length(args)) {
  eval(parse(text = args[i]))
}

dataset <- GetDataReducedAA("SELECT * FROM binary_labeling_activity")

mc_itr = 1000
con <- get_con("sfp_AcpS")

CrossValidation <- function(enzyme, num_leave_out, alpha.0.prior.coef, alpha.1.prior.coef, p1, to.plot) {
  alpha.0.prior <- SetPriorReducedAA(alpha.0.prior.coef, NUM_CLASS)
  alpha.1.prior <- SetPriorReducedAA(alpha.1.prior.coef, NUM_CLASS)
  dataset_X <- dataset$feature[dataset$data[, enzyme] != -1, ]
  dataset_Y <- dataset$data[dataset$data[, enzyme] != -1, enzyme]
  cv_result <- leave_n_out_cv(dataset_X, dataset_Y, alpha.1.prior, alpha.0.prior, p1, mc_itr, num_leave_out)
  roc.xy <- GetRocXy(cv_result$prob, cv_result$Y)
  if (to.plot == TRUE) {
    dbWriteTable(con, value = roc.xy , name = paste0("ROC_", enzyme), overwrite = TRUE, append = FALSE)
  } else {
    auc <- GetAuc(roc.xy$x, roc.xy$y)
    query <- sprintf("INSERT INTO cv_find_params VALUES ('%s', %f, %f, %f, %f)", enzyme, alpha.0.prior.coef, alpha.1.prior.coef, p1, auc)
    data <- dbGetQuery(conn=con, statement=query)
  }
}

if (TOPLOT == 1) {
  CrossValidation('sfp', 3, 10, 0.1, 1e-5, TRUE)
  CrossValidation('AcpS', 2, 10, 1, 1e-5, TRUE)
  CrossValidation('PfAcpH', 3, 10, 10, 0.5, TRUE)
} else {
  CrossValidation('sfp', 3, ALPHA0, ALPHA1, P1, FALSE)
  CrossValidation('AcpS', 2, ALPHA0, ALPHA1, P1, FALSE)
  CrossValidation('PfAcpH', 3, ALPHA0, ALPHA1, P1, FALSE)
}

dbDisconnect(con)
