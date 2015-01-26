root_path <- "../.."
datafile <- "../../data/all_TS_data.csv"

org.data <- read.csv(datafile, stringsAsFactors=F)

#### plot cumulative number of hits v.s. number of rounds
pdf("cumulative_number_of_hits_vs_number_of_rounds.pdf", width=10, height=10)
TS <- c(1,3,4,5)
par(mfrow=c(2,2))
# sfp specific labeling
counts <- rep(0, 4)
for (i in 1:dim(org.data)[1]) {
    if (org.data[i, 'sfp'] == 1 && org.data[i, 'AcpS'] == 0) {
        counts[which(TS == org.data[i, 'TS'])] <- counts[which(TS == org.data[i, 'TS'])] + 1
    }
}
cumu_counts = rep(0,4)
for (i in 1:length(counts)) {
    cumu_counts[i] = sum(counts[1:i])
}
barplot(cumu_counts, main="sfp specific labeling cumulative", xlab="number of rounds")
# sfp specific labeling and unlabeling
counts <- rep(0, 4)
for (i in 1:dim(org.data)[1]) {
    if (org.data[i, 'type1'] == 1) {
        counts[which(TS == org.data[i, 'TS'])] <- counts[which(TS == org.data[i, 'TS'])] + 1
    }
}
cumu_counts = rep(0,4)
for (i in 1:length(counts)) {
    cumu_counts[i] = sum(counts[1:i])
}
barplot(cumu_counts, main="sfp specific labeling and unlabeling cumulative", xlab="number of rounds")
# AcpS specific labeling
counts <- rep(0, 4)
for (i in 1:dim(org.data)[1]) {
    if (org.data[i, 'AcpS'] == 1 && org.data[i, 'sfp'] == 0) {
        counts[which(TS == org.data[i, 'TS'])] <- counts[which(TS == org.data[i, 'TS'])] + 1
    }
}
cumu_counts = rep(0,4)
for (i in 1:length(counts)) {
    cumu_counts[i] = sum(counts[1:i])
}
barplot(cumu_counts, main="AcpS specific labeling cumulative", xlab="number of rounds")
# AcpS specific labeling and unlabeling
counts <- rep(0, 4)
for (i in 1:dim(org.data)[1]) {
    if (org.data[i, 'type2'] == 1) {
        counts[which(TS == org.data[i, 'TS'])] <- counts[which(TS == org.data[i, 'TS'])] + 1
    }
}
cumu_counts = rep(0,4)
for (i in 1:length(counts)) {
    cumu_counts[i] = sum(counts[1:i])
}
barplot(cumu_counts, main="AcpS specific labeling and unlabeling cumulative", xlab="number of rounds")
dev.off()

#### plot histogram of length of type 1 hits
pdf("histogram_of_length_type1.pdf", width=10, height=10)
lengths <- c()
for (i in 1:dim(org.data)[1]) {
    if (org.data[i, 'type1'] == 1) {
        lengths <- c(lengths, nchar(org.data[i, 'sequence']))
    }
}
hist(lengths, breaks=10)
dev.off()

#### plot diverse hits, this code produces data and to generate pdf files, run python script
datafile <- paste(root_path, "/data/TS0.csv", sep='')
source(paste(root_path,"/src/NB_Greedy_library/NB_utility.R",sep=''))
source(paste(root_path,"/src/NB_Greedy_library/NB_interface.R",sep=''))
classfile <- paste(root_path,"/data/Reduced_AA_Alphabet.csv",sep='')
num_rec <- 50
add_ins <- 5
gamma_0 <- 100
gamma_1 <- 1
prior <- 1e-4
nL <- 19
nR <- 19
maxL <- 8
maxR <- 8
minL <- 8
minR <- 8
S.Pos <- nL
data_org <- read.csv(datafile, stringsAsFactors=F)
AAclass <- read.csv(classfile, stringsAsFactors=F)
train_data <- getFeatures(data_org, AAclass, nL, nR)
X <- train_data[, 1:(nL+nR)]
Y <- train_data[, 'sfp']
rec_pool <- generate_recommendation_MAP_old(X, Y, AAclass, S.Pos, num_rec, maxL, maxR, minL, minR, gamma_0, gamma_1, prior, add_ins, 1)
pool_class <- convert_feature_to_class(rec_pool$rec, AAclass, nL, nR, maxL, maxR)
fout <- file("pool.txt")
writeLines(pool_class, fout)
close(fout)
# mutation
rec_mutate <- mutate_rec(data_org[,'nterm'], data_org[,'cterm'], AAclass, num_rec, 2)
mutate_seqs <- convert_to_seqs(rec_mutate$nterm, rec_mutate$cterm, maxL, maxR)
mutate_class <- convert_to_class(mutate_seqs, AAclass)
fout <- file("mutate.txt")
writeLines(mutate_class, fout)
close(fout)
# naive
rec_naive <- naive_rec(X, Y, AAclass, S.Pos, nL, nR, gamma_0, gamma_1, prior, 100, 1e4, num_rec, maxL, maxR)
naive_seqs <- convert_to_seqs(rec_naive$nterm, rec_naive$cterm, maxL, maxR)
naive_class <- convert_to_class(naive_seqs, AAclass)
fout <- file("naive.txt")
writeLines(naive_class, fout)
close(fout)

#### plot performance benchmark of POOL, mutation and naive method
datafile <- paste(root_path, "/data/TS0.csv", sep='')
all_datafile <- paste(root_path, "/data/all_TS_data.csv", sep='')
source(paste(root_path,"/src/NB_Greedy_library/NB_utility.R",sep=''))
source(paste(root_path,"/src/NB_Greedy_library/NB_interface.R",sep=''))
source("utility_for_plot.R")
classfile <- paste(root_path,"/data/Reduced_AA_Alphabet.csv",sep='')
AAclass <- read.csv(classfile, stringsAsFactors=F)
num_rec <- 50
add_ins <- 1
gamma_0 <- 100
gamma_1 <- 0.05
prior <- 1e-4
nL <- 19
nR <- 19
maxL <- 12
maxR <- 12
minL <- 10
minR <- 10
itr <- 500
S.Pos <- nL
all_data_org <- read.csv(all_datafile, stringsAsFactors=F)[1:8]
all_train_data <- getFeatures(all_data_org, AAclass, nL, nR)
X_sfp <- all_train_data[all_train_data[,'sfp'] != -1, 1:(nL+nR)]
Y_sfp <- all_train_data[all_train_data[,'sfp'] != -1, 'sfp']
data_org <- read.csv(datafile, stringsAsFactors=F)
train_data <- getFeatures(data_org, AAclass, nL, nR)
X <- train_data[, 1:(nL+nR)]
Y <- train_data[, 'sfp']
X <- rbind(X, X_sfp[1:268,])
Y <- c(Y, Y_sfp[1:268])
# pool
rec_pool <- generate_recommendation_MAP_old(X, Y, AAclass, S.Pos, num_rec, maxL, maxR, minL, minR, gamma_0, gamma_1, prior, add_ins, 1)
X_pool <- rec_pool$rec
# mutation
rec_mutate <- mutate_rec(data_org[,'nterm'], data_org[,'cterm'], AAclass, num_rec, 2)
X_mutate <- convert_term_to_feature(rec_mutate$nterm, rec_mutate$cterm, AAclass, nL, nR)
# naive
rec_naive <- naive_rec(X, Y, AAclass, S.Pos, nL, nR, gamma_0, gamma_1, prior, 100, 1e4, num_rec, maxL, maxR)
X_naive <- convert_term_to_feature(rec_naive$nterm, rec_naive$cterm, AAclass, nL, nR)

prob_pool <- Naive_Bayes(X_sfp, Y_sfp, X_pool, AAclass, S.Pos, nL, nR, gamma_0, gamma_1, prior, itr)
prob_mutate <- Naive_Bayes(X_sfp, Y_sfp, X_mutate, AAclass, S.Pos, nL, nR, gamma_0, gamma_1, prior, itr)
prob_naive <- Naive_Bayes(X_sfp, Y_sfp, X_naive, AAclass, S.Pos, nL, nR, gamma_0, gamma_1, prior, itr)
# plot
plot_data <- data.frame(x=1:length(prob_pool), POOL=get_P(prob_pool), Mutation=get_P(prob_mutate), Naive=get_P(prob_naive))
require(ggplot2)
require(reshape2)
dd <- melt(plot_data, id=c("x"))
g <- ggplot(dd) + geom_line(aes(x=x, y=value, colour=variable))
g <- g + labs(x="Number of peptides to sample", y="P(at least one peptide is a hit)", title="POOL works better than alternative methods")
g <- g + scale_color_discrete(name="Method")
ggsave("benchmark.pdf", g)

