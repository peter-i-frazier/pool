rm(list = ls())
source('norm_util.R')
root.path <- "../.."
org.data <- read.csv(sprintf("%s/data/all_TS_raw_reading.csv", root.path), stringsAsFactors=F)
result.path <- sprintf("%s/src/normalization/norm_result", root.path)
treatment <- "sfp_1"
# preprocess
reduced.data <- org.data[!is.na(org.data[, treatment]),]
j.class <- unique(reduced.data[, 'TS'])
j.class <- j.class[order(j.class)]
for (ts in j.class) {
    min <- min(reduced.data[reduced.data[, 'TS'] == ts, treatment])
    max <- max(reduced.data[reduced.data[, 'TS'] == ts, treatment])
    reduced.data[reduced.data[, 'TS'] == ts, treatment] <- (reduced.data[reduced.data[, 'TS'] == ts, treatment] - min) / (max - min)
}
unique.seqs <- unique(reduced.data[, 'seq'])
dup.data <- c()
for (seq in unique.seqs) {
    sub.data <- reduced.data[reduced.data[, 'seq'] == seq, c('TS', 'spot', 'seq', treatment)]
    if (nrow(sub.data) > 1) {
        dup.data <- rbind(dup.data, sub.data)
    }
}
# end preprocess
normalize(reduced.data, dup.data, treatment, result.path=result.path)

# plot
data <- read.csv(sprintf('%s/%s_norm_table.csv', result.path, treatment), stringsAsFactors=F)
dup.data <- read.csv(sprintf('%s/%s_dup_table.csv', result.path, treatment), stringsAsFactors=F)
# plot density of theta
pdf(sprintf("%s/%s_theta_hist.pdf", result.path, treatment))
legend.name <- c('TS1','TS2', 'TS3', 'TS4', 'TS5')
plot(c(-0.1,1.2), c(0,15), type='n', xlab='theta', ylab='density')
co <- c('blue', 'red', 'green', 'purple', 'cyan')
for (i in j.class) {
    # d <- density(data[data[,'TS'] == i, 'theta'])
    # lines(d, col=co[i])
    d <- hist(data[data[,'TS'] == i, 'theta'], breaks=30, plot=F)
    lines(d$mids, d$density, col=co[i])
}
legend('topright', legend=legend.name[j.class], col=co[j.class], lwd=1)
dev.off()
# plot density of epsilon for duplicates
pdf(sprintf("%s/%s_eps_hist.pdf", result.path, treatment))
hist(dup.data[,'eps'], breaks=50, xlab='epsilon', ylab='count', xlim=c(-0.5, 0.5))
dev.off()
