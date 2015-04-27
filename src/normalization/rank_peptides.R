rm(list=ls())

root_path = "../.."
org.data <- read.csv(paste(root_path, "/data/all_TS_raw_reading.csv", sep=''), stringsAsFactors=F)

rank_by_measure <- function(table.list) {
    mu.vec <- rep(0, 4)
    sd.vec <- rep(0, 4)
    sfp.1 <- table.list[[1]]
    sfp.2 <- table.list[[2]]
    AcpS.1 <- table.list[[3]]
    AcpS.2 <- table.list[[4]]
    mu.vec[1] <- mean(sfp.1[, 'theta'])
    sd.vec[1] <- sd(sfp.1[, 'theta'])
    #mu.vec[2] <- mean(sfp.1[, 'theta'] - sfp.2[, 'theta'])
    mu.vec[2] <- mean(sfp.2[, 'theta'])
    #sd.vec[2] <- sd(sfp.1[, 'theta'] - sfp.2[, 'theta'])
    sd.vec[2] <- sd(sfp.2[, 'theta'])
    mu.vec[3] <- mean(AcpS.1[, 'theta'])
    sd.vec[3] <- sd(AcpS.1[, 'theta'])
    mu.vec[4] <- mean(AcpS.1[, 'theta'] - AcpS.2[, 'theta'])
    sd.vec[4] <- sd(AcpS.1[, 'theta'] - AcpS.2[, 'theta'])
# measure =1: non-specific sfp-type hit
    measure.1 <- sfp.1[, 1:5]
    colnames(measure.1)[5] <- 'measure'
    for (n in 1:nrow(measure.1)) {
        #v1 <- org.data[org.data[, 'seq'] == measure.1[n, 'seq'], 'sfp_2'][1]
        v1 <- sfp.1[sfp.1[, 'seq'] == measure.1[n, 'seq'], 'theta'][1]
        #v2 <- sfp.2[sfp.2[, 'seq'] == measure.1[n, 'seq'], 'theta'][1]
        #measure.1[n, 'measure'] <- (v1 - mu.vec[1]) / sd.vec[1] + 2 * (v1 - v2 - mu.vec[2]) / sd.vec[2]
        #measure.1[n, 'measure'] <- (v1 - mu.vec[1]) / sd.vec[1] - 0.9 * (v2 - mu.vec[2]) / sd.vec[2]
        measure.1[n, 'measure'] <- v1
    }
    measure.1 <- measure.1[order(measure.1[, 'measure'], decreasing=T),]
#         =2: non-specific AcpS-type hit
    measure.2 <- AcpS.1[, 1:5]
    colnames(measure.2)[5] <- 'measure'
    for (n in 1:nrow(measure.2)) {
        v3 <- AcpS.1[AcpS.1[, 'seq'] == measure.2[n, 'seq'], 'theta'][1]
        v4 <- AcpS.2[AcpS.2[, 'seq'] == measure.2[n, 'seq'], 'theta'][1]
        measure.2[n, 'measure'] <- (v3 - mu.vec[3]) / sd.vec[3] + (v3 - v4 - mu.vec[4]) / sd.vec[4]
    }
    measure.2 <- measure.2[order(measure.2[, 'measure'], decreasing=T),]
#         =3: specific sfp-type hit
    measure.3 <- AcpS.1[, 1:5]
    colnames(measure.3)[5] <- 'measure'
    for (n in 1:nrow(measure.3)) {
        v1 <- sfp.1[sfp.1[, 'seq'] == measure.3[n, 'seq'], 'theta'][1]
        v2 <- sfp.2[sfp.2[, 'seq'] == measure.3[n, 'seq'], 'theta'][1]
        v3 <- AcpS.1[AcpS.1[, 'seq'] == measure.3[n, 'seq'], 'theta'][1]
        measure.3[n, 'measure'] <- (v1 - mu.vec[1]) / sd.vec[1]
                                 + (v1 - v2 - mu.vec[2]) / sd.vec[2]
                                 - (v3 - mu.vec[3]) / sd.vec[3]
    }
    measure.3 <- measure.3[order(measure.3[, 'measure'], decreasing=T),]
#         =4: specific AcpS-type hit
    measure.4 <- AcpS.1[, 1:5]
    colnames(measure.4)[5] <- 'measure'
    for (n in 1:nrow(measure.4)) {
        v1 <- sfp.1[sfp.1[, 'seq'] == measure.4[n, 'seq'], 'theta'][1]
        v3 <- AcpS.1[AcpS.1[, 'seq'] == measure.4[n, 'seq'], 'theta'][1]
        v4 <- AcpS.2[AcpS.2[, 'seq'] == measure.4[n, 'seq'], 'theta'][1]
        measure.4[n, 'measure'] <- (v3 - mu.vec[3]) / sd.vec[3]
                                 + (v3 - v4 - mu.vec[4]) / sd.vec[4]
                                 - (v1 - mu.vec[1]) / sd.vec[1]
    }
    measure.4 <- measure.4[order(measure.4[, 'measure'], decreasing=T),]
    return(list(measure.1, measure.2, measure.3, measure.4))
}

slice.table <- function(org.table, unit.length, dest.folder, filename) {
    org.table <- org.table[org.table[, 'TS'] != 2,]
    for (i in 0:floor(nrow(org.table) / unit.length))
        write.csv(org.table[(i * unit.length + 1):((i + 1) * unit.length),], sprintf("%s/%s_%d_to_%d.csv", dest.folder, filename, i * unit.length + 1, (i+1) * unit.length), row.names=F)
    write.csv(org.table[(floor(nrow(org.table) / unit.length) * unit.length + 1):nrow(org.table),], sprintf("%s/%s_%d_to_%d.csv", dest.folder, filename, floor(nrow(org.table) / unit.length) * unit.length + 1, nrow(org.table)), row.names=F)
}

##### script
measure.list <- rank_by_measure(list(read.csv('sfp_1_norm_result_table.csv', stringsAsFactors=F),
                                     read.csv('sfp_2_norm_result_table.csv', stringsAsFactors=F),
                                     read.csv('AcpS_1_norm_result_table.csv', stringsAsFactors=F),
                                     read.csv('AcpS_2_norm_result_table.csv', stringsAsFactors=F)))
filenames <- c('measure1', 'measure2', 'measure3', 'measure4')
for (i in 1:4) {
    slice.table(measure.list[[i]], 100, 'rank_table', filenames[i])
    write.csv(measure.list[[i]], sprintf('rank_table/%s.csv', filenames[i]))
}
