# Author: Jialei Wang
# Generate data for benchmark plot

root <- "src/"
source(paste0(root, "core/dependency.R"))
ResolveDependency(root)

args <- commandArgs(TRUE)
num_recom <- as.numeric(args[1])
label_name <- args[2]
notlabel_name <- args[3]
specific_name <- args[4]

data <- c()
if (specific_name == 'sfp_specific') {
  end_idx <- 2
} else {
  end_idx <- 3
}
for (i in 0:end_idx) {
  data <- rbind(data, read.csv(sprintf('data/TS%d.csv', i)))
}
AA_lookup <- read.csv('data/reduced_AA_lookup.csv')
AA_lookup[,'AA'] = as.character(AA_lookup[,'AA'])
train_dataset <- GetDataReducedAA(data, AA_lookup)

# pool generate_recommendation
prior.df <- data.frame(sfp = c(10, 0.1, 1e-4), AcpS = c(10, 1, 1e-4), sfp_specific = c(10, 1, 1e-4), acps_specific = c(10, 1, 1e-4))
label.alpha.0 <- SetPriorReducedAA(prior.df[1, label_name], NUM_CLASS)
label.alpha.1 <- SetPriorReducedAA(prior.df[2, label_name], NUM_CLASS)
label.p1 <- prior.df[3, label_name]

notlabel.alpha.0 <- SetPriorReducedAA(prior.df[1, notlabel_name], NUM_CLASS)
notlabel.alpha.1 <- SetPriorReducedAA(prior.df[2, notlabel_name], NUM_CLASS)
notlabel.p1 <- prior.df[3, notlabel_name]

specific_alpha_0 <- SetPriorReducedAA(prior.df[1, specific_name], NUM_CLASS)
specific_alpha_1 <- SetPriorReducedAA(prior.df[2, specific_name], NUM_CLASS)
specific_p1 <- prior.df[3, specific_name]

label_idx <- train_dataset$data[, label_name] != -1
notlabel_idx <- train_dataset$data[, notlabel_name] != -1

pool.recoms <- GenRecomSetMAPNewOrthogonal(train_dataset$feature[label_idx, ], train_dataset$data[label_idx, label_name],
train_dataset$feature[notlabel_idx, ], train_dataset$data[notlabel_idx, notlabel_name],
label.alpha.1, label.alpha.0, label.p1,
notlabel.alpha.1, notlabel.alpha.0, notlabel.p1,
1000, num_recom, minL, maxL, minR, maxR)
# end pool

mutate.recoms <- GenRecomSetMutate(train_dataset$feature[(train_dataset$data[, label_name] == 1) & (train_dataset$data[, notlabel_name] == 0), ], 4, 1, num_recom)

predict.optimize.recoms <- GenRecomSetPredictOptimizeOrthogonal(train_dataset$feature[label_idx, ], train_dataset$data[label_idx, label_name],
train_dataset$feature[notlabel_idx, ], train_dataset$data[notlabel_idx, notlabel_name],
label.alpha.1, label.alpha.0, label.p1,
notlabel.alpha.1, notlabel.alpha.0, notlabel.p1,
1000, num_recom, minL, maxL, minR, maxR)

all_data <- c()
for (i in 0:5) {
  all_data <- rbind(all_data, read.csv(sprintf('data/TS%d.csv', i)))
}
true_dataset <- GetDataReducedAA(all_data, AA_lookup)
true_label_idx <- true_dataset$data[, label_name] != -1
true_notlabel_idx <- true_dataset$data[, notlabel_name] != -1
true_specific_idx <- true_dataset$data[, specific_name] != -1


pool.quality <- ComputeProbImprovOfSet(pool.recoms, true_dataset$feature[true_specific_idx,], true_dataset$data[true_specific_idx, specific_name], specific_alpha_1, specific_alpha_0, specific_p1, 1000)
mutate.quality <- ComputeProbImprovOfSet(mutate.recoms, true_dataset$feature[true_specific_idx,], true_dataset$data[true_specific_idx, specific_name], specific_alpha_1, specific_alpha_0, specific_p1, 1000)
predict.optimize.quality <- ComputeProbImprovOfSet(predict.optimize.recoms, true_dataset$feature[true_specific_idx,], true_dataset$data[true_specific_idx, specific_name], specific_alpha_1, specific_alpha_0, specific_p1, 1000)

quality.result <- data.frame(pool = pool.quality, mutate = mutate.quality, predict_optimize = predict.optimize.quality)
write.csv(quality.result, paste0(specific_name, '_simulation_quality.csv'), row.names=FALSE)
# con <- get_con("sfp_AcpS")
# dbWriteTable(con, value = quality.result, name = paste0("benchmark_", type_name), overwrite = TRUE, append = FALSE, row.names=FALSE)
#
pool.nterm <- c()
pool.cterm <- c()
mutate.nterm <- c()
mutate.cterm <- c()
predict.optimize.nterm <- c()
predict.optimize.cterm <- c()
for (i in 1:num_recom) {
  pool.list <- GenRandomSeqFromFeature(pool.recoms[i,], AA_lookup)
  pool.nterm <- c(pool.nterm, pool.list$nterm)
  pool.cterm <- c(pool.cterm, pool.list$cterm)
  mutate.list <- GenRandomSeqFromFeature(mutate.recoms[i,], AA_lookup)
  mutate.nterm <- c(mutate.nterm, mutate.list$nterm)
  mutate.cterm <- c(mutate.cterm, mutate.list$cterm)
  predict.optimize.list <- GenRandomSeqFromFeature(predict.optimize.recoms[i,], AA_lookup)
  predict.optimize.nterm <- c(predict.optimize.nterm, predict.optimize.list$nterm)
  predict.optimize.cterm <- c(predict.optimize.cterm, predict.optimize.list$cterm)
}
seq.table <- data.frame(pool_nterm = pool.nterm, pool_cterm = pool.cterm, mutation_nterm = mutate.nterm, mutation_cterm = mutate.cterm,
naive_nterm = predict.optimize.nterm, naive_cterm = predict.optimize.cterm, stringsAsFactors=FALSE)
write.csv(seq.table, paste0(specific_name, '_simulation_seq.csv'), row.names=FALSE)
# dbWriteTable(con, value = seq.table, name = paste0('benchmark_seq_', type_name), overwrite = TRUE, append = FALSE, row.names=FALSE)
#
# original.seq.table <- cbind(train_dataset$data[,c('nterm', 'cterm')], 1*((train_dataset$data[, label_name] == 1) & (train_dataset$data[, notlabel_name] == 0)))
# colnames(original.seq.table) <- c('original_nterm', 'original_cterm', 'hit')
# dbWriteTable(con, value = original.seq.table, name = paste0('original_seq_', type_name), overwrite = TRUE, append = FALSE, row.names=FALSE)
#
# dbDisconnect(con)
