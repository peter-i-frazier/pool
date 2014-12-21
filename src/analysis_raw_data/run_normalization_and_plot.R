rm(list=ls())
cutoff <- 2.0
data_file <- '/Users/jialeiwang/Downloads/raw_data.csv'
DATA <- read.csv(data_file)
source('analysis.R')

norm_data <- c()
for (i in 1:4) {
    norm_data <- cbind(norm_data, normalize_data(DATA[,i], cutoff))
}
write.csv(norm_data, 'norm_data.csv')

par(mfrow=c(2,2))
for (i in 1:4) {
    hist(norm_data[,i], breaks=70)
}


