rm(list=ls())

root_path <- "/fs/home/jw865/peptide-catalysis"

# parameters that can change
num_recom <- 105

type1_gamma_0 <- 100
type1_gamma_1 <- 0.01
type1_prior <- 1e-5

type2_gamma_0 <- 10
type2_gamma_1 <- 1
type2_prior <- 1e-5

add_ins <- 1
nL <- 19
nR <- 19
max_left <- 9
max_right <- 9
min_left <- 4
min_right <- 4
serine_position <- nL
source(paste(root_path,"/src/NB_Greedy_library/Naive_Bayes_util.R",sep=''))
source(paste(root_path,"/src/NB_Greedy_library/utility.R",sep=''))
source(paste(root_path,"/src/NB_Greedy_library/Opt_Search_util.R",sep=''))
AAclass <- data.frame(read.csv(paste(root_path, "/data/Reduced_AA_Alphabet.csv", sep=''), header = T, as.is = T, sep = ','))
data.org <- read.csv(paste(root_path, "/data/2014_06_03_orthogonal_labeling_data/Training_Set_Cumulative.csv", sep=''),as.is=T)
data <- getFeatures(data.org, AAclass, nL, nR)
trainX_type1 <- data[,1:(nL+nR)]
trainY_type1 <- data[,'type1_hit']
trainX_type2 <- trainX_type1
trainY_type2 <- data[,'type2_hit']

nAA <- length(unique(as.numeric(AAclass)))
# find type 1 hit
    which_type <- 'type 1 hit old method'
    X <- trainX_type1
    Y <- trainY_type1
    gamma_0 <- type1_gamma_0
    gamma_1 <- type1_gamma_1
    prior <- type1_prior
# find type 2 hit
#    which_type <- 'type 2 hit old method'
#    X <- trainX_type2
#    Y <- trainY_type2
#    gamma_0 <- type1_gamma_0
#    gamma_1 <- type1_gamma_1
#    prior <- type1_prior

recom_list <- MaxProbSearchMAPOld(X, Y, AAclass, serine_position, num_recom, max_left, max_right, min_left, min_right, gamma_0, gamma_1, add_ins)

# validate prob of recom_list
#Parallel Computation Packages
require(snow)
require(doParallel)
#=================================================================================
#Parallel Computation Setting
# hosts <- rep('whale',40)
# hosts <- c(rep("whale",52),rep("ahab",24),rep("flask",24),rep("starbuck",24))
# hosts <- c(rep("whale",30),rep("ahab",30),rep("flask",30),rep("starbuck",30))
hosts <- c(rep("whale",10),rep("ahab",10),rep("flask",10),rep("starbuck",10))
# hosts <- c(rep("ishmael",30),rep("stubb",30),rep("daggoo",30),rep("tashtego",30))
# hosts <- c(rep("ishmael",10),rep("stubb",10),rep("daggoo",10),rep("tashtego",10))
cl <- makeCluster(hosts, type = "SOCK")
registerDoParallel(cl)

iteration <- 500
alpha <- Dirichlet_Parameter(X, Y, AAclass, gamma_0, gamma_1)

predict_mat <- foreach (n = 1:iteration, .init = c(), .combine = rbind, .inorder=F, .errorhandling='pass') %dopar% {
    theta <- getTheta_MC(alpha = alpha, classlist = AAclass)
    NB_predict(recom_list, theta, serine_position, max_left, max_right, prior)
}
recom_prob <- colMeans(predict_mat)
print (which_type)
print (add_ins)

stopCluster(cl)

