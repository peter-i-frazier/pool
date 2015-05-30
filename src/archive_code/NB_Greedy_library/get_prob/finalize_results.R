root_path <- "/fs/home/jw865/peptide-catalysis"
datafile <- paste(root_path,"/data/all_TS_data.csv",sep='')
size_table <- 2573
first_N <- size_table

#import libraries
source(paste(root_path,"/src/NB_Greedy_library/NB_utility.R",sep=''))
source(paste(root_path,"/src/NB_Greedy_library/NB_interface.R",sep=''))
data_org <- read.csv(datafile, header = T, as.is = T)
prob_table <- c()
for (i in 1:size_table) {
    one_result <- t(read.csv(paste(root_path, "/temp/get_prob/prob_", i, ".csv", sep=''))[,2])
    prob_table <- rbind(prob_table, one_result)
}
colnames(prob_table) <- c('sfp', 'AcpS', 'PfAcpH', 'org_idx')
# calculate prob of combinations
type_1_prob <- prob_table[, 'sfp'] * (1 - prob_table[, 'AcpS']) * prob_table[, 'PfAcpH']
type_2_prob <- prob_table[, 'AcpS'] * (1 - prob_table[, 'sfp']) * prob_table[, 'PfAcpH']
both_label_prob <- prob_table[, 'sfp'] * prob_table[, 'AcpS']
neither_label_prob <- (1 - prob_table[, 'sfp']) * (1 - prob_table[, 'AcpS'])
# rank them
type_1_table <- cbind(data_org[prob_table[order(type_1_prob, decreasing=T)[1:first_N], 'org_idx'],], type_1_prob[order(type_1_prob, decreasing=T)[1:first_N]])
colnames(type_1_table) <- c('sfp', 'AcpS', 'PfAcpH', 'type1', 'type2', 'nterm', 'cterm', 'seq', 'TS', 'spot', 'type_1_prob')
type_2_table <- cbind(data_org[prob_table[order(type_2_prob, decreasing=T)[1:first_N], 'org_idx'],], type_2_prob[order(type_2_prob, decreasing=T)[1:first_N]])
colnames(type_2_table) <- c('sfp', 'AcpS', 'PfAcpH', 'type1', 'type2', 'nterm', 'cterm', 'seq', 'TS', 'spot', 'type_2_prob')
both_label_table <- cbind(data_org[prob_table[order(both_label_prob, decreasing=T)[1:first_N], 'org_idx'],], both_label_prob[order(both_label_prob, decreasing=T)[1:first_N]])
colnames(both_label_table) <- c('sfp', 'AcpS', 'PfAcpH', 'type1', 'type2', 'nterm', 'cterm', 'seq', 'TS', 'spot', 'both_label_prob')
prob_table_out <- cbind(data_org[prob_table[, 'org_idx'],], type_1_prob, type_2_prob, both_label_prob)
colnames(prob_table_out) <- c('sfp', 'AcpS', 'PfAcpH', 'type1', 'type2', 'nterm', 'cterm', 'seq', 'TS', 'spot', 'type_1_prob', 'type_2_prob', 'both_label_prob')
# pick top 10
neither_label_top <- cbind(data_org[prob_table[order(neither_label_prob, decreasing=T)[1:10], 'org_idx'],], neither_label_prob[order(neither_label_prob, decreasing=T)[1:10]])
colnames(neither_label_top) <- c('sfp', 'AcpS', 'PfAcpH', 'type1', 'type2', 'nterm', 'cterm', 'seq', 'TS', 'spot', 'neither_label_prob')
type_1_top <- c()
for (ts in c(1,3,4,5)) {
    count <- 0
    idx <- 0
    while (count <= 10) {
        idx <- idx + 1
        if (type_1_table[idx, 'TS'] == ts) {
            type_1_top <- rbind(type_1_top, type_1_table[idx,])
            count <- count + 1
        }
    }
}
colnames(type_1_top) <- colnames(type_1_table)
type_2_top <- c()
for (ts in c(1,3,4,5)) {
    count <- 0
    idx <- 0
    while (count <= 10) {
        idx <- idx + 1
        if (type_2_table[idx, 'TS'] == ts) {
            type_2_top <- rbind(type_2_top, type_2_table[idx,])
            count <- count + 1
        }
    }
}
colnames(type_2_top) <- colnames(type_2_table)
both_label_top <- c()
for (ts in c(1,3,4,5)) {
    count <- 0
    idx <- 0
    while (count <= 10) {
        idx <- idx + 1
        if (both_label_table[idx, 'TS'] == ts) {
            both_label_top <- rbind(both_label_top, both_label_table[idx,])
            count <- count + 1
        }
    }
}
colnames(both_label_top) <- colnames(both_label_table)
# write to file
write.csv(type_1_top, 'type_1_recommendation.csv')
write.csv(type_2_top, 'type_2_recommendation.csv')
write.csv(both_label_top, 'both_label_recommendation.csv')
write.csv(neither_label_top, 'neither_label_recommendation.csv')
write.csv(prob_table_out, 'prob_lookup_table.csv')


