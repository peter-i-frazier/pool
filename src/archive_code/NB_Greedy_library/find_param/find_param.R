# This file needs root_path, datafile, outcome_name, para_idx as input in command line
# Code block below is only for testing, must disable when submit to nbs!
# root_path = "/fs/home/jw865/peptide-catalysis"
# datafile <- paste(root_path,"/data/whole_experiment_data.csv",sep='')
# outcome_name <- "sfp"
# para_idx <- 300
# num_fold <- 3
# Code block end

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

gamma_0_list <- c(1000, 500, 100, 50, 10, 5, 1)
gamma_1_list <- c(0.01, 0.05, 0.1, 0.5, 1, 10)
prior_list <- c(1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.2, 0.3, 0.4, 0.5)
parameters <- group_parameters(gamma_0_list, gamma_1_list, prior_list)
colnames(parameters) <- c('gamma0','gamma1','prior')

gamma_0 <- parameters[para_idx, 'gamma0']
gamma_1 <- parameters[para_idx, 'gamma1']
prior <- parameters[para_idx, 'prior']
nL <- 19
nR <- 19
S.Pos <- 19
maxL <- 19
maxR <- 19
itr <- 500

#import libraries
source(paste(root_path,"/src/NB_Greedy_library/NB_utility.R",sep=''))
source(paste(root_path,"/src/NB_Greedy_library/NB_interface.R",sep=''))

#Specify path/to/classlist file here
classfile <- paste(root_path,"/data/Reduced_AA_Alphabet.csv",sep='')

data_org <- read.csv(datafile, header = T, as.is = T)
AAclass <- read.csv(classfile, header = T, as.is = T)
train_data <- getFeatures(data_org, AAclass, nL, nR)

original_X <- train_data[,1:(nL+nR)]
original_Y <- train_data[,outcome_name]
X <- original_X[original_Y != -1,]
Y <- original_Y[original_Y != -1]
cv_result <- fold_cv(X, Y, AAclass, S.Pos, nL, nR, gamma_0, gamma_1, prior, itr, num_fold) 
auc <- cv_result$auc
# cv_result <- loocv(X, Y, AAclass, S.Pos, nL, nR, gamma_0, gamma_1, prior, itr) 
# xy <- cv_result$axis
# area <- cv_result$area
# 
# plot(x = xy$x, y = xy$y, type = 'l', xlim = c(0,1), ylim = c(0,1), xlab = "false positive rate", ylab = "true positive rate")
# title(main = )
# dev.off()

write(auc, 'auc')
print ("auc")
print (auc)
print ("gamma 0")
print (gamma_0)
print ("gamma 1")
print (gamma_1)
print ("prior")
print (prior)
print (outcome_name)
