rm(list=ls())
source("normalization_util.R")

rank_seqs_by_measure <- function(data, theta_table, letter_number_table, measure, coef_1, coef_2, coef_3) {
# measure =1: non-specific sfp-type hit
#         =2: non-specific AcpS-type hit
#         =3: specific sfp-type hit
#         =4: specific AcpS-type hit
    seq_table <- cbind(theta_table[, 'unique_seqs'], rep(NA, dim(theta_table)[1]))
    colnames(seq_table) <- c('seq', 'measure')
# compute measure
    for (n in 1:dim(theta_table)[1]) {
        if (measure == 1) {
            seq_table[n, 'measure'] <- coef_1 * theta_table[n, 'sfp_1'] + coef_2 * (theta_table[n, 'sfp_1'] - theta_table[n, 'sfp_2'])
        } else if (measure == 2) {
            seq_table[n, 'measure'] <- coef_1 * theta_table[n, 'AcpS_1'] + coef_2 * (theta_table[n, 'AcpS_1'] - theta_table[n, 'AcpS_2'])
        } else if (measure == 3) {
            seq_table[n, 'measure'] <- coef_1 * theta_table[n, 'sfp_1'] + coef_2 * (theta_table[n, 'sfp_1'] - theta_table[n, 'sfp_2']) - coef_3 * theta_table[n, 'AcpS_1']
        } else if (measure == 4) {
            seq_table[n, 'measure'] <- coef_1 * theta_table[n, 'AcpS_1'] + coef_2 * (theta_table[n, 'AcpS_1'] - theta_table[n, 'AcpS_2']) - coef_3 * theta_table[n, 'sfp_1']
        }
    }
# rank according to measure
    seq_table <- seq_table[order(seq_table[, 'measure'], decreasing=T),]
# construct lookup table
    lookup_table <- c()
    for (n in 1:dim(seq_table)[1]) {
        part_table <- data[data[, 'seq'] == seq_table[n, 'seq'], c('seq', 'TS', 'spot')]
        num_replica <- dim(part_table)[1]
        part_table <- cbind(part_table, rep(NA, num_replica), rep(NA, num_replica), rep(seq_table[n, 'measure'], num_replica))
        colnames(part_table) <- c('seq', 'TS', 'spot', 'row', 'col', 'measure')
        for (i in 1:num_replica) {
            row_letter <- strsplit(part_table[i, 'spot'], split='')[[1]][1]
            part_table[i, 'col'] <- as.numeric(paste(strsplit(part_table[i, 'spot'], split='')[[1]][-1], collapse='')) - 1
            part_table[i, 'row'] <- letter_number_table[1, row_letter]
        }
        lookup_table <- rbind(lookup_table, part_table)
    }
    return (lookup_table)
}

root_path = "../.."
org.data <- read.csv(paste(root_path, "/data/all_TS_raw_reading.csv", sep=''), stringsAsFactors=F)
theta_table <- read.csv("theta_table.csv", stringsAsFactor=F)
letter_number_table <- read.csv("letter_number_lookup.csv", stringsAsFactors=F)
std_sfp <- sd(theta_table[, 'sfp_1'])
std_sfp_diff <- sd(theta_table[, 'sfp_1'] - theta_table[, 'sfp_2'])
# theta_table that excludes TS2
truncated_theta_table <- theta_table[!is.na(theta_table[, 'AcpS_1']),]
std_AcpS <- sd(truncated_theta_table[, 'AcpS_1'])
std_AcpS_diff <- sd(truncated_theta_table[, 'AcpS_1'] - truncated_theta_table[, 'AcpS_2'])
write.csv(rank_seqs_by_measure(org.data, theta_table, letter_number_table, 1,  1/std_sfp, 1/std_sfp_diff, NA),
          "rank_table/rank_by_measure_1.csv", row.names=F)
write.csv(rank_seqs_by_measure(org.data, truncated_theta_table, letter_number_table, 2,  1/std_AcpS, 1/std_AcpS_diff, NA),
          "rank_table/rank_by_measure_2.csv", row.names=F)
write.csv(rank_seqs_by_measure(org.data, truncated_theta_table, letter_number_table, 3,  1/std_sfp, 1/std_sfp_diff, 1/std_AcpS),
          "rank_table/rank_by_measure_3.csv", row.names=F)
write.csv(rank_seqs_by_measure(org.data, truncated_theta_table, letter_number_table, 4,  1/std_AcpS, 1/std_AcpS_diff, 1/std_sfp),
          "rank_table/rank_by_measure_4.csv", row.names=F)
