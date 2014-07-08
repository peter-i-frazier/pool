# This file needs root_path, datafile, outcome_name, para_idx as input in command line
# Code block below is only for testing, must disable when submit to nbs!
# root_path = "/fs/home/jw865/peptide-catalysis"
# datafile <- paste(root_path,"/data/2014_06_03_orthogonal_labeling_data/Training_Set_Cumulative.csv",sep='')
# outcome_name <- "PfAcpH"
# para_idx <- 1
# Code block end
load(paste(root_path, "/data/all_data_for_training.RData", sep=''))

X <- trainX_sfp
Y <- trainY_sfp

group_parameters <- function(v1, v2, v3) {
    l1 <- length(v1)
    l2 <- length(v2)
    l3 <- length(v3)
    group <- matrix(-1,nrow=l1*l2*l3,ncol=3)
    ct <- 0
    for (i in 1:l1)
        for (j in 1:l2)
            for (k in 1:l3) {
                ct <- ct + 1
                group[ct,1] <- v1[i]
                group[ct,2] <- v2[j]
                group[ct,3] <- v3[k]
            }
    return (group)
}

args <- commandArgs(TRUE)
for(i in 1:length(args))
{
	eval(parse(text = args[i]))
}

gamma_0_list <- c(1000, 500, 100, 50, 10, 1)
gamma_1_list <- c(0.01, 0.05, 0.1, 0.5, 1)
prior_list <- c(1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.3, 0.5)
parameters <- group_parameters(gamma_0_list, gamma_1_list, prior_list)
colnames(parameters) <- c('gamma0','gamma1','prior')

Gamma_0 <- parameters[para_idx, 'gamma0']
Gamma_1 <- parameters[para_idx, 'gamma1']
prior.positive <- parameters[para_idx, 'prior']
filename <- 'loocv'
nL <- 19
nR <- 19
S.Pos <- 19
maxL <- 19
maxR <- 19
predIter <- 500

#import libraries
source(paste(root_path,"/src/NB_Greedy_library/Naive_Bayes_util.R",sep=''))
source(paste(root_path,"/src/NB_Greedy_library/utility.R",sep=''))
source(paste(root_path,"/src/NB_Greedy_library/Opt_Search_util.R",sep=''))
#Specify path/to/classlist file here
classfile <- paste(root_path,"/data/Reduced_AA_Alphabet.csv",sep='')

# data_org <- data.frame(read.csv(datafile, header = T, as.is = T, sep = ','))
AAclass <- data.frame(read.csv(classfile, header = T, as.is = T, sep = ','))
# train_data <- getFeatures(data_org, AAclass, nL, nR)
# training part has 2 cases: if we train PfAcpH, it's the prob give the peptide is labeled, this training set is a subset of original training set
# X <- train_data[,1:(nL+nR)]
# N <- dim(X)[1]
# if (outcome_name == 'PfAcpH') {
#     truth.table <- (((train_data[,'sfp_specific']==1) + (train_data[,'AcpS_specific']==1)) >= 1)
#     train.X <- train_data[truth.table,1:(nL+nR)]
#     train.Y <- train_data[truth.table,'PfAcpH']
#     sub.N <- length(train.Y)
#     print("length of train.Y")
#     print(sub.N)
#     sub.prob <- c()
#     for (n in 1:sub.N) {
#         sub.prob <- c(sub.prob, Naive_Bayes(train.X[-n,], train.Y[-n], train.X[n,], AAclass, S.Pos, maxL, maxR, Gamma_0, Gamma_1, prior.positive, predIter))
#     }
#     ROC_plot(sub.prob, train.Y, filename)
#     xy <- ROC_xy_output(sub.prob, train.Y)
#     prob <- Naive_Bayes(train.X, train.Y, X, AAclass, S.Pos, maxL, maxR, Gamma_0, Gamma_1, prior.positive, predIter)
#     ct <- 0
#     for (i in 1:N) {
#         if (truth.table[i]) {
#             ct <- ct + 1
#             prob[i] <- sub.prob[ct]
#         }
#     }
# } else {
#     Y <- train_data[,outcome_name]
#     print(length(Y))
#     prob <- c()
#     for (n in 1:N) {
#         prob <- c(prob, Naive_Bayes(X[-n,], Y[-n], X[n,], AAclass, S.Pos, maxL, maxR, Gamma_0, Gamma_1, prior.positive, predIter))
#     }
#     ROC_plot(prob, Y, filename)
#     xy <- ROC_xy_output(prob, Y)
# }

prob <- c()
for (n in 1:N) {
   prob <- c(prob, Naive_Bayes(X[-n,], Y[-n], X[n,], AAclass, S.Pos, maxL, maxR, Gamma_0, Gamma_1, prior.positive, predIter))
}
ROC_plot(prob, Y, filename)
xy <- ROC_xy_output(prob, Y)

auc <- AUC(rev(xy$x), rev(xy$y))
write(auc, 'auc')
write.csv(prob, 'prob.csv')
print ("auc")
print (auc)
print ("gamma 0")
print (Gamma_0)
print ("gamma 1")
print (Gamma_1)
print ("prior")
print (prior.positive)
print (outcome_name)
