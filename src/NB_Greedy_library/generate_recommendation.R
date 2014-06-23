rm(list=ls())
root_path <- "/Users/jialeiwang/peptide-catalysis"
data_path <- "/data/2014_06_03_orthogonal_labeling_data/Training_Set_Cumulative.csv"
source(paste(root_path,"/src/NB_Greedy_library/Naive_Bayes_util.R",sep=''))
source(paste(root_path,"/src/NB_Greedy_library/utility.R",sep=''))
source(paste(root_path,"/src/NB_Greedy_library/Opt_Search_util.R",sep=''))
data_org <- data.frame(read.csv(paste(root_path, data_path, sep=''), header = T, as.is = T, sep =','))
AAclass <- data.frame(read.csv(paste(root_path, "/data/Reduced_AA_Alphabet.csv", sep=''), header = T, as.is = T, sep = ','))
nL <- 19
nR <- 19
maxL <- 10
maxR <- 10
minL <- 4
minR <- 4
Nlib <- 10
iteration <- 10
S.Pos <- nL
nF <- (nL + nR)
nAA <- length(unique(as.numeric(AAclass)))
train_data <- getFeatures(data_org, AAclass, nL, nR)

### old method
# old_outcome <- 'type1_hit'
# no_recom <- 500
# Gamma_0 <- 100
# Gamma_1 <- 0.01
# prior <- 1e-4
# trainX <- train_data[,1:nF]
# trainY <- train_data[,old_outcome]
# recom_list <- c()
# for (pep_idx in 1:no_recom) {
#     print (pep_idx)
#     print (old_outcome)
#     pep_lib <- gen_peptide_lib(Nlib, nAA, nF, maxL, maxR, minL, minR, S.Pos, colnames(trainX))
#     prob <- Naive_Bayes(trainX, trainY, pep_lib, AAclass, S.Pos, maxL, maxR, Gamma_0, Gamma_1, prior, iteration)
#     new_pep <- select_new_recom(recom_list, pep_lib, prob, pep_idx)
#     recom_list <- rbind(recom_list, new_pep)
#     # add to trainX
#     trainX <- rbind(trainX, new_pep)
#     trainY <- c(trainY, 0)
# }

### new method
no_recom <- 500
truth.table <- (((train_data[,'sfp_specific']==1) + (train_data[,'AcpS_specific']==1)) >= 1)

new_step1_outcome <- 'sfp_specific'
step1_Gamma_0 <- 100
step1_Gamma_1 <- 0.01
step1_prior <- 1e-4
step1_X <- train_data[,1:nF]
step1_Y <- train_data[,new_step1_outcome]

new_step2_outcome <- 'PfAcpH'
step2_Gamma_0 <- 100
step2_Gamma_1 <- 0.01
step2_prior <- 1e-4
step2_X <- train_data[truth.table,1:nF]
step2_Y <- train_data[truth.table,new_step2_outcome]

recom_list <- c()
for (pep_idx in 1:no_recom) {
    print (pep_idx)
    print (new_step1_outcome)
    pep_lib <- gen_peptide_lib(Nlib, nAA, nF, maxL, maxR, minL, minR, S.Pos, colnames(step1_X))
    step1_prob <- Naive_Bayes(step1_X, step1_Y, pep_lib, AAclass, S.Pos, maxL, maxR, step1_Gamma_0, step1_Gamma_1, step1_prior, iteration)
    step2_prob <- Naive_Bayes(step2_X, step2_Y, pep_lib, AAclass, S.Pos, maxL, maxR, step2_Gamma_0, step2_Gamma_1, step2_prior, iteration)
    prob <- step1_prob * step2_prob
    new_pep <- select_new_recom(recom_list, pep_lib, prob, pep_idx)
    recom_list <- rbind(recom_list, new_pep)
    # add to step1_X, step2_X
    step1_X <- rbind(step1_X, new_pep)
    step1_Y <- c(step1_Y,0)
    step2_X <- rbind(step2_X, new_pep)
    step2_Y <- c(step2_Y,0)
}
