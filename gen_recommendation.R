# Author: Jialei Wang
# Generate recommendation

root <- "src/"
source(paste0(root, "core/dependency.R"))
ResolveDependency(root)

args <- commandArgs(TRUE)
round <- as.numeric(args[1])
num_recom <- as.numeric(args[2])
label_name <- args[3]
notlabel_name <- args[4]
unlabel_name <- args[5]
fname <- args[6]
print(sprintf("Compose peptides for %s", fname))

data <- c()
for (i in 0:as.numeric(round)) {
  data <- rbind(data, read.csv(sprintf('data/TS%d.csv', i)))
}
AA_lookup <- read.csv('data/reduced_AA_lookup.csv')
AA_lookup[,'AA'] = as.character(AA_lookup[,'AA'])

train_dataset <- GetDataReducedAA(data, AA_lookup)
# pool generate_recommendation
prior.df <- data.frame(sfp = c(10, 0.1, 1e-4), AcpS = c(10, 1, 1e-4), PfAcpH = c(10, 10, 0.5))
label.alpha.0 <- SetPriorReducedAA(prior.df[1, label_name], NUM_CLASS)
label.alpha.1 <- SetPriorReducedAA(prior.df[2, label_name], NUM_CLASS)
label.p1 <- prior.df[3, label_name]
label_idx <- train_dataset$data[, label_name] != -1

if (notlabel_name != 'null') {
  notlabel.alpha.0 <- SetPriorReducedAA(prior.df[1, notlabel_name], NUM_CLASS)
  notlabel.alpha.1 <- SetPriorReducedAA(prior.df[2, notlabel_name], NUM_CLASS)
  notlabel.p1 <- prior.df[3, notlabel_name]
  notlabel_idx <- train_dataset$data[, notlabel_name] != -1
} else {
  notlabel.alpha.0 <- c()
  notlabel.alpha.1 <- c()
  notlabel.p1 <- c()
  notlabel_idx <- c()
}

unlabel.alpha.0 <- SetPriorReducedAA(prior.df[1, unlabel_name], NUM_CLASS)
unlabel.alpha.1 <- SetPriorReducedAA(prior.df[2, unlabel_name], NUM_CLASS)
unlabel.p1 <- prior.df[3, unlabel_name]
unlabel_idx <- train_dataset$data[, unlabel_name] != -1

pool.recoms <- GenRecomSetNew(train_dataset$feature[label_idx, ], train_dataset$data[label_idx, label_name],
train_dataset$feature[notlabel_idx, ], train_dataset$data[notlabel_idx, notlabel_name],
train_dataset$feature[unlabel_idx, ], train_dataset$data[unlabel_idx, unlabel_name],
label.alpha.1, label.alpha.0, label.p1,
notlabel.alpha.1, notlabel.alpha.0, notlabel.p1,
unlabel.alpha.1, unlabel.alpha.0, unlabel.p1,
1000, num_recom, minL, maxL, minR, maxR)
## end pool
nterms <- c()
cterms <- c()
seqs <- c()
for (i in 1:num_recom) {
  peptide <- (GenRandomSeqFromFeature(pool.recoms[i,], AA_lookup))
  nterms <- c(nterms, peptide$nterm[1])
  cterms <- c(cterms, peptide$cterm[1])
  seqs <- c(seqs, paste0(peptide$nterm[1], 'S', peptide$cterm[1]))
}
df <- data.frame(nterm=nterms, cterm=cterms, seq=seqs)
print(df)
write.csv(df, fname, row.names=FALSE)
