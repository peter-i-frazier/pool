#Procedures for analyzing data

#0. Read in Data
rm(list=ls())
#DATA <- read.csv('../../data/2014_08_08_orthogonal_labeling/sfp_vs_AcpS.csv')
DATA <- read.csv('~/Desktop/data.csv')
source('analysis.R')

#1. Normalize Raw Data:
#        Input: raw data y
#        1. Find median of raw data M
#        2. Find 25% quantile of raw data L
#        3. Find 25% quantile of standard normal z
#        4. std of raw data \sigma = (M-L) / -z
#        5. normalized data y' =  (y-M) / (\sigma)
#        Output: y'
#
#2. Specify threshold by observing plot of normalized data
#h <- hist(normalize_data(DATA[,2]), 100, xaxp=c(-5,15,50))
#cursor_x <- rep(4, 2)
#cursor_y <- c(0,15)
#lines(cursor_x, cursor_y, col="red")

#3. For each peptide x, it has four indicator values:
#        x_1 = 1 if labeled by sfp
#        x_2 = 0 if unlabeled by PfAcpH (in sfp membrane)
#        x_3 = 1 if labeled by Acps
#        x_4 = 0 if unlabeled by PfAcpH (in Acps membrane)
#        Let x = (x_1, x_2, x_3, x_4)
#threshold <- c(1670, 1310, 1420, 1470)
#threshold <- c(1824, 1425, 1557, 1624)
#indicator <- peptide_indicator_raw(DATA, threshold)

#4. x can be put into one of the 7 categories:
#        1. not labeled:                             (0, 0, 0, 0)
#        2. labeled by both, unlabeled:              (1, 0, 1, 0)
#        3. labeled by both, fail to unlabel:        (1, 1, 1, 1) 
#        4. only labeled by sfp, unlabeled:          (1, 0, 0, 0)
#        5. only labeled by sfp, fail to unlabel:    (1, 1, 0, 0)
#        6. only labeled by Acps, unlabeled:         (0, 0, 1, 0)
#        7. only labeled by Acps, fail to unlabel:   (0, 0, 1, 1)
#

# script start
normal_data <- c()
for (i in 1:2) {
    normal_data <- cbind(normal_data, normalize_data(DATA[,i],2.0))
}
write.csv(normal_data,'normal.csv')
#M <- c(1470.38, 1078.21, 1010.46, 1064.98)
#for (i in 1:4) {
#    normal_data <- cbind(normal_data, normalize_data_manual(DATA[,i],M[i]))
#}

#####################
## Find hit by categorizing peptides
##We are interested in 4 & 6
#threshold <- c(1.5, 3, 1.5, 3)
#indicator <- peptide_indicator(normal_data, threshold)
#cat <- categorize_peptide(indicator)
#sfp_code_idx <- c(which(cat==4), which(cat==5))
#sfp_spot_idx <- convert_to_spot_idx(sfp_code_idx)
#sfp_category <- cbind(sfp_spot_idx, 1-indicator[sfp_code_idx,2], normal_data[sfp_code_idx,1:4])
#colnames(sfp_category) <- c('Index', 'unlabeling', 'sfp normalized value', 'sfp-PfAcpH normalized value',  'AcpS normalized value', 'AcpS-PfAcpH normalized value')
#write.csv(sfp_category, 'sfp_category.csv')
#
#AcpS_code_idx <- c(which(cat==6), which(cat==7))
#AcpS_spot_idx <- convert_to_spot_idx(AcpS_code_idx)
#AcpS_category <- cbind(AcpS_spot_idx, 1-indicator[AcpS_code_idx,4], normal_data[AcpS_code_idx,1:4])
#colnames(AcpS_category) <- c('Index', 'unlabeling', 'sfp normalized value', 'sfp-PfAcpH normalized value', 'AcpS normalized value', 'AcpS-PfAcpH normalized value')
#write.csv(AcpS_category, 'AcpS_category.csv')
#####################
#
#
#####################
## Find hit by taking diff
#labeling_cutoff <- rep(1.5, 2)
#delta_cutoff <- rep(0.5, 2)
##M <- c(1493, 1134.4, 1192, 1207.7)
##sigma <- c(190, 120, 125, 125)
##normal_data <- normalize_data_manual(DATA, M, sigma)
#
#sfp_labeling <- 1 * (normal_data[,1] > labeling_cutoff[1])
#AcpS_labeling <- 1 * (normal_data[,3] > labeling_cutoff[2])
#both_labeling <- sfp_labeling * AcpS_labeling
#only_sfp_labeling <- sfp_labeling - both_labeling
#only_AcpS_labeling <- AcpS_labeling - both_labeling
#delta_sfp <- normal_data[,1] - normal_data[,2]
#delta_AcpS <- normal_data[,3] - normal_data[,4]
#sfp_unlabel <- 1 * (delta_sfp > delta_cutoff[1])
#AcpS_unlabel <- 1 * (delta_AcpS > delta_cutoff[2])
#
#sfp_code_idx <- which(only_sfp_labeling==1)
#sfp_spot_idx <- convert_to_spot_idx(sfp_code_idx)
#sfp_diff <- cbind(sfp_spot_idx, sfp_unlabel[sfp_code_idx], normal_data[sfp_code_idx,1:4], delta_sfp[sfp_code_idx])
#colnames(sfp_diff) <- c('Index', 'unlabeling', 'sfp normalized value', 'sfp-PfAcpH normalized value',  'AcpS normalized value', 'AcpS-PfAcpH normalized value', 'sfp diff')
#write.csv(sfp_diff, 'sfp_diff.csv')
#
#AcpS_code_idx <- which(only_AcpS_labeling==1)
#AcpS_spot_idx <- convert_to_spot_idx(AcpS_code_idx)
#AcpS_diff <- cbind(AcpS_spot_idx, AcpS_unlabel[AcpS_code_idx], normal_data[AcpS_code_idx,1:4], delta_AcpS[AcpS_code_idx])
#colnames(AcpS_diff) <- c('Index', 'unlabeling', 'sfp normalized value', 'sfp-PfAcpH normalized value', 'AcpS normalized value', 'AcpS-PfAcpH normalized value', 'AcpS diff')
#write.csv(AcpS_diff, 'AcpS_diff.csv')
#####################
#
#
#####################
## Analyze Lori's result
#sfp_data <- read.csv('Lori_result_sfp.csv')
#AcpS_data <- read.csv('Lori_result_AcpS.csv')
#sfp_spot_idx <- as.character(sfp_data[,1])
#sfp_code_idx <- convert_to_code_idx(sfp_spot_idx)
#sfp_analysis <- cbind(sfp_data, normal_data[sfp_code_idx,1:4])
#colnames(sfp_analysis) <- c('Index', 'unlabeling', 'sfp normalized value', 'sfp-PfAcpH normalized value',  'AcpS normalized value', 'AcpS-PfAcpH normalized value')
#write.csv(sfp_analysis, 'Lori_result_sfp_analysis.csv')
#
#AcpS_spot_idx <- as.character(AcpS_data[,1])
#AcpS_code_idx <- convert_to_code_idx(AcpS_spot_idx)
#AcpS_analysis <- cbind(AcpS_data, normal_data[AcpS_code_idx,1:4])
#colnames(AcpS_analysis) <- c('Index', 'unlabeling', 'sfp normalized value', 'sfp-PfAcpH normalized value', 'AcpS normalized value', 'AcpS-PfAcpH normalized value')
#write.csv(AcpS_analysis, 'Lori_result_AcpS_analysis.csv')
#####################





