rm(list=ls())

load('RDATA_result10/sfp_addin_10.RData')

pred_prob <- Naive_Bayes(train_data[,1:nF], train_data[,old_outcome], recom_list, AAclass, S.Pos, maxL, maxR, Gamma_0, Gamma_1, prior, iteration)

source('utility.R')

writePep(recom_list, S.Pos, pred_prob, AAclass, 'type1_addin_10')
