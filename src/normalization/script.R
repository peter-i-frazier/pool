rm(list=ls())
require(mclust)
data_file <- 'all_TS_raw_reading.csv'
data <- read.csv(data_file, stringsAsFactors=F)
pep_seq = 'DSLEFIASKLA'
#pep_seq = 'DALEFIASKLA'
#pep_seq = 'GDSLSWLLRLLN'

# pdf plot
par(mfrow=c(2,3))
for (which_TS in 1:4) {
    v <- data[data[,'TS']==which_TS, 'sfp.step.1']
    mc <- Mclust(v)
    plot(mc, what='density', xlab='x', ylab='density')
    lines(density(v), col='red')
    subdata <- data[data[,'TS'] == which_TS,]
    points(subdata[subdata[,'seq'] == pep_seq, 'sfp.step.1'], 0, col='red')
    value = subdata[subdata[,'seq'] == pep_seq, 'sfp.step.1']
    print ('norm')
    dens <- densityMclust(v)
    print (cdfMclust(dens, value)$y)
    #print (summary(mc, parameters=T))
}


# # cdf plot
# par(mfrow=c(2,3))
# for (which_TS in 1:5) {
#     v <- data[data[,'TS']==which_TS, 'sfp.step.1']
#     mc <- Mclust(v)
#     subdata <- data[data[,'TS'] == which_TS,]
#     value = subdata[subdata[,'seq'] == pep_seq, 'sfp.step.1']
#     dens <- densityMclust(v)
#     plot(ecdf(v), main='CDF', xlab='x', ylab='P(X <= x)', col='red')
#     lines(cdfMclust(dens))
#     points(value, cdfMclust(dens, value)$y, col='red', pch=19)
# }
