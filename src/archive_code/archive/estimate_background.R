rm(list=ls())
root_path = "../.."
org.data <- read.csv(paste(root_path, "/data/all_TS_raw_reading.csv", sep=''), stringsAsFactors=F)
background.data <- read.csv('background.csv', stringsAsFactors=F)
for (r in 1:nrow(background.data))
    background.data[r, 'val'] <- org.data[org.data[,'TS'] == background.data[r, 'TS'] & 
                                          org.data[,'spot'] == background.data[r, 'spot'], 
                                          background.data[r, 'treatment']]
mean.table <- c()
for (ts in 1:5)
    for (k in c('sfp_1', 'sfp_2', 'AcpS_1', 'AcpS_2'))
        if (!(ts==2 && (k=='AcpS_1' || k=='AcpS_2')))
            mean.table <- rbind(mean.table, data.frame(TS=ts, treatment=k, 
                                     mean=mean(background.data[background.data[,'TS']==ts & background.data[,'treatment']==k, 'val']), 
                                     stringsAsFactors=F))
for (r in 1:nrow(background.data))
    background.data[r, 'diff'] <- (background.data[r, 'val'] - mean.table[mean.table[, 'TS'] == background.data[r, 'TS'] & mean.table[, 'treatment'] == background.data[r, 'treatment'], 'mean']) / background.data[r, 'val']
