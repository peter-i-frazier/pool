# Author: Jialei Wang
# Created: 07.22.2015
# Generate data for benchmark plot

rm(list = ls())
root <- "/fs/home/jw865/remote_deployment/ucsd_reversible_labeling/"
source(paste0(root, "core/dependency.R"))
ResolveDependency(root)

#args <- commandArgs(TRUE)
#label_name <- args[1]
#notlabel_name <- args[2]
#type_name <- args[3]

label_name <- 'sfp'
notlabel_name <- 'AcpS'
type_name <- 'type1'
num_recom <- 100

train_dataset <- GetDataReducedAA("SELECT * FROM binary_labeling_activity WHERE TS <= 2")

# pool generate_recommendation
prior.df <- data.frame(sfp = c(10, 0.1, 1e-4), AcpS = c(10, 1, 1e-4), PfAcpH = c(10, 10, 0.5))
label.alpha.0 <- SetPriorReducedAA(prior.df[1, label_name], NUM_CLASS)
label.alpha.1 <- SetPriorReducedAA(prior.df[2, label_name], NUM_CLASS)
label.p1 <- prior.df[3, label_name]

notlabel.alpha.0 <- SetPriorReducedAA(prior.df[1, notlabel_name], NUM_CLASS)
notlabel.alpha.1 <- SetPriorReducedAA(prior.df[2, notlabel_name], NUM_CLASS)
notlabel.p1 <- prior.df[3, notlabel_name]

PfAcpH.alpha.0 <- SetPriorReducedAA(prior.df[1, 'PfAcpH'], NUM_CLASS)
PfAcpH.alpha.1 <- SetPriorReducedAA(prior.df[2, 'PfAcpH'], NUM_CLASS)
PfAcpH.p1 <- prior.df[3, 'PfAcpH']

label_idx <- train_dataset$data[, label_name] != -1
notlabel_idx <- train_dataset$data[, notlabel_name] != -1
PfAcpH_idx <- train_dataset$data[, 'PfAcpH'] != -1

pool.recoms <- GenRecomSetMAPNew(train_dataset$feature[label_idx, ], train_dataset$data[label_idx, label_name],
                                 train_dataset$feature[notlabel_idx, ], train_dataset$data[notlabel_idx, notlabel_name],
                                 train_dataset$feature[PfAcpH_idx, ], train_dataset$data[PfAcpH_idx, 'PfAcpH'],
                                 label.alpha.1, label.alpha.0, label.p1,
                                 notlabel.alpha.1, notlabel.alpha.0, notlabel.p1,
                                 PfAcpH.alpha.1, PfAcpH.alpha.0, PfAcpH.p1,
                                 1000, num_recom, minL, maxL, minR, maxR)
# end pool

mutate.recoms <- GenRecomSetMutate(train_dataset$feature[train_dataset$data[, type_name] == 1, ], 4, 1, num_recom)

predict.optimize.recoms <- GenRecomSetPredictOptimize(train_dataset$feature[label_idx, ], train_dataset$data[label_idx, label_name],
                                                      train_dataset$feature[notlabel_idx, ], train_dataset$data[notlabel_idx, notlabel_name],
                                                      train_dataset$feature[PfAcpH_idx, ], train_dataset$data[PfAcpH_idx, 'PfAcpH'],
                                                      label.alpha.1, label.alpha.0, label.p1,
                                                      notlabel.alpha.1, notlabel.alpha.0, notlabel.p1,
                                                      PfAcpH.alpha.1, PfAcpH.alpha.0, PfAcpH.p1,
                                                      1000, num_recom, minL, maxL, minR, maxR)

# calculate quality of recoms
GetQuality <- function(label.X, label.Y, not.label.X, not.label.Y, unlabel.X, unlabel.Y,
                       label.alpha.1.prior, label.alpha.0.prior, label.p1,
                       not.label.alpha.1.prior, not.label.alpha.0.prior, not.label.p1,
                       unlabel.alpha.1.prior, unlabel.alpha.0.prior, unlabel.p1,
                       test.X, num.mc) {
  label.prob <- Predict(label.X, label.Y, test.X, label.alpha.1.prior, label.alpha.0.prior, label.p1, num.mc)
  not.label.prob <- Predict(not.label.X, not.label.Y, test.X, not.label.alpha.1.prior, not.label.alpha.0.prior, not.label.p1, num.mc)
  unlabel.prob <- Predict(unlabel.X, unlabel.Y, test.X, unlabel.alpha.1.prior, unlabel.alpha.0.prior, unlabel.p1, num.mc)
  prob <- label.prob * (1 - not.label.prob) * unlabel.prob

  num_rec <- length(prob)
  quality <- rep(0, num_rec)
  prod <- 1
  for (i in 1:num_rec) {
      prod <- prod * (1 - prob[i])
      quality[i] <- 1 - prod
  }
  return (quality)
}

true_dataset <- GetDataReducedAA("SELECT * FROM binary_labeling_activity")
true_label_idx <- true_dataset$data[, label_name] != -1
true_notlabel_idx <- true_dataset$data[, notlabel_name] != -1
true_PfAcpH_idx <- true_dataset$data[, 'PfAcpH'] != -1

pool.quality <- GetQuality(true_dataset$feature[true_label_idx, ], true_dataset$data[true_label_idx, label_name],
                           true_dataset$feature[true_notlabel_idx, ], true_dataset$data[true_notlabel_idx, notlabel_name],
                           true_dataset$feature[true_PfAcpH_idx, ], true_dataset$data[true_PfAcpH_idx, 'PfAcpH'],
                           label.alpha.1, label.alpha.0, label.p1,
                           notlabel.alpha.1, notlabel.alpha.0, notlabel.p1,
                           PfAcpH.alpha.1, PfAcpH.alpha.0, PfAcpH.p1,
                           pool.recoms, 1000)

mutate.quality <- GetQuality(true_dataset$feature[true_label_idx, ], true_dataset$data[true_label_idx, label_name],
                           true_dataset$feature[true_notlabel_idx, ], true_dataset$data[true_notlabel_idx, notlabel_name],
                           true_dataset$feature[true_PfAcpH_idx, ], true_dataset$data[true_PfAcpH_idx, 'PfAcpH'],
                           label.alpha.1, label.alpha.0, label.p1,
                           notlabel.alpha.1, notlabel.alpha.0, notlabel.p1,
                           PfAcpH.alpha.1, PfAcpH.alpha.0, PfAcpH.p1,
                           mutate.recoms, 1000)

predict.optimize.quality <- GetQuality(true_dataset$feature[true_label_idx, ], true_dataset$data[true_label_idx, label_name],
                           true_dataset$feature[true_notlabel_idx, ], true_dataset$data[true_notlabel_idx, notlabel_name],
                           true_dataset$feature[true_PfAcpH_idx, ], true_dataset$data[true_PfAcpH_idx, 'PfAcpH'],
                           label.alpha.1, label.alpha.0, label.p1,
                           notlabel.alpha.1, notlabel.alpha.0, notlabel.p1,
                           PfAcpH.alpha.1, PfAcpH.alpha.0, PfAcpH.p1,
                           predict.optimize.recoms, 1000)

quality.result <- data.frame(pool = pool.quality, mutate = mutate.quality, predict_optimize = predict.optimize.quality)
con <- get_con("sfp_AcpS")
dbWriteTable(con, value = quality.result, name = paste0("benchmark_", type_name), overwrite = TRUE, append = FALSE, row.names=FALSE)

lookup.table <- QuerySql("sfp_AcpS", "SELECT * FROM reduced_AA")
pool.nterm <- c()
pool.cterm <- c()
mutate.nterm <- c()
mutate.cterm <- c()
predict.optimize.nterm <- c()
predict.optimize.cterm <- c()
for (i in 1:num_recom) {
  pool.list <- GenRandomSeqFromFeature(pool.recoms[i,], lookup.table)
  pool.nterm <- c(pool.nterm, pool.list$nterm)
  pool.cterm <- c(pool.cterm, pool.list$cterm)
  mutate.list <- GenRandomSeqFromFeature(mutate.recoms[i,], lookup.table)
  mutate.nterm <- c(mutate.nterm, mutate.list$nterm)
  mutate.cterm <- c(mutate.cterm, mutate.list$cterm)
  predict.optimize.list <- GenRandomSeqFromFeature(predict.optimize.recoms[i,], lookup.table)
  predict.optimize.nterm <- c(predict.optimize.nterm, predict.optimize.list$nterm)
  predict.optimize.cterm <- c(predict.optimize.cterm, predict.optimize.list$cterm)
}
seq.table <- data.frame(pool_nterm = pool.nterm, pool_cterm = pool.cterm, mutation_nterm = mutate.nterm, mutation_cterm = mutate.cterm,
naive_nterm = predict.optimize.nterm, naive_cterm = predict.optimize.cterm, stringsAsFactors=FALSE)
dbWriteTable(con, value = seq.table, name = 'benchmark_seq', overwrite = TRUE, append = FALSE, row.names=FALSE)

original.seq.table <- train_dataset$data[,c('nterm', 'cterm')]
colnames(original.seq.table) <- c('original_nterm', 'original_cterm')
dbWriteTable(con, value = original.seq.table, name = 'original_seq', overwrite = TRUE, append = FALSE, row.names=FALSE)

dbDisconnect(con)
