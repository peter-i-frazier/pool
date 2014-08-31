rm(list=ls())
#Parallel Computation Packages
require(snow)
require(doParallel)
#=================================================================================
#Parallel Computation Setting
# hosts <- rep('whale',40)
# hosts <- c(rep("whale",52),rep("ahab",24),rep("flask",24),rep("starbuck",24))
# hosts <- c(rep("whale",30),rep("ahab",30),rep("flask",30),rep("starbuck",30))
hosts <- c(rep("ishmael",30),rep("stubb",30),rep("daggoo",30),rep("tashtego",30))
cl <- makeCluster(hosts, type = "SOCK")
registerDoParallel(cl)

root_path = "/fs/home/jw865/peptide-catalysis"
#==============================for training new method
# load(paste(root_path, "/data/all_data_for_training.RData", sep=''))

gamma_0_list <- c(500, 100, 50, 10, 1)
gamma_1_list <- c(0.01, 0.05, 0.1, 0.5, 1)
prior_list <- c(1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.5)
# gamma_0_list <- c(100)
# gamma_1_list <- c(0.5)
# prior_list <- c(1e-4)

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
AAclass <- read.csv(classfile, header = T, as.is = T, sep = ',')
#========================== for training type1 & type2 old method
data.org <- data.frame(read.csv(paste(root_path, "/data/2014_06_03_orthogonal_labeling_data/Training_Set_Cumulative.csv", sep=''),header=T, as.is=T, sep=","))
data <- getFeatures(data.org, AAclass, nL, nR)
trainX_type1 <- data[,1:(nL+nR)]
trainY_type1 <- data[,'type1_hit']
trainX_type2 <- trainX_type1
trainY_type2 <- data[,'type2_hit']
#=============================================
outcome_name <- "type1_hit"
X <- trainX_type1
Y <- trainY_type1
# X <- trainX_type2
# Y <- trainY_type2


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
            prob <- foreach (n = 1:dim(X)[1], .init = c(), .combine = c, .inorder = T, .errorhandling = 'pass') %dopar% {
                    Naive_Bayes(X[-n,], Y[-n], X[n,], AAclass, S.Pos, maxL, maxR, gamma_0_list[i], gamma_1_list[j], prior_list[k], predIter)
            }
            xy <- ROC_xy_output(prob, Y)
            auc <- c(auc, AUC(rev(xy$x), rev(xy$y)))
        }
    }
}

max.index <- which.max(auc)
print (index[max.index,])
print (auc[max.index])
