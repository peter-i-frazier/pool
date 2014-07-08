rm(list=ls())
#Parallel Computation Packages
require(multicore)
require(snow)
require(doParallel)
require(MCMCpack)
require(Rlab)
#=================================================================================
#Parallel Computation Setting
# hosts <- rep('whale',40)
# hosts <- c(rep("whale",52),rep("ahab",24),rep("flask",24),rep("starbuck",24))
hosts <- c(rep("whale",30),rep("ahab",30),rep("flask",30),rep("starbuck",30))
# hosts <- c(rep("ishmael",30),rep("stubb",30),rep("daggoo",30),rep("tashtego",30))
cl <- makeCluster(hosts, type = "SOCK")
registerDoParallel(cl)

root_path = "/fs/home/jw865/peptide-catalysis"
load(paste(root_path, "/data/all_data_for_training.RData", sep=''))

outcome_name <- "sfp"
X <- trainX_sfp
Y <- trainY_sfp

gamma_0_list <- c(1000, 500, 100, 50, 10, 1)
gamma_1_list <- c(0.01, 0.05, 0.1, 0.5, 1)
prior_list <- c(1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.3, 0.5)

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
AAclass <- data.frame(read.csv(classfile, header = T, as.is = T, sep = ','))

auc <- c()
index <- matrix(-1, nrow=length(gamma_0_list)*length(gamma_1_list)*length(prior_list), ncol=3)
count <- 0
for (i in 1:length(gamma_0_list)) {
    for (j in 1:length(gamma_1_list)) {
        for (k in 1:length(prior_list)) {
            count <- count + 1
            print (count)
            print (outcome_name)
            index[count,] <- c(i,j,k)
            prob <- foreach (n = 1:dim(X)[1], .init = c(), .combine = c, .packages = 'Rlab', .inorder = T, .errorhandling = 'pass') %dopar% {
                    Naive_Bayes(X[-n,], Y[-n], X[n,], AAclass, S.Pos, maxL, maxR, gamma_0_list[i], gamma_1_list[j], prior_list[k], predIter)
            }
            xy <- ROC_xy_output(prob, Y)
            auc <- c(auc, AUC(rev(xy$x), rev(xy$y)))
        }
    }
}
