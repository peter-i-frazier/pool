rm(list=ls())

root_path <- "/fs/home/jw865/peptide-catalysis"

# parameters that can change
no_recom <- 60

type1_gamma_0 <- 100
type1_gamma_1 <- 0.01
type1_prior <- 1e-5

type2_gamma_0 <- 10
type2_gamma_1 <- 1
type2_prior <- 1e-5

add_ins <- 1
nL <- 19
nR <- 19
maxL <- 9
maxR <- 9
minL <- 4
minR <- 4
Nlib <- 1e7
S.Pos <- nL
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
#    which_type <- 'type 1 hit'
#    X <- trainX_type1
#    Y <- trainY_type1
#    gamma_0 <- type1_gamma_0
#    gamma_1 <- type1_gamma_1
#    prior <- type1_prior
# find type 2 hit
    which_type <- 'type 2 hit'
    X <- trainX_type2
    Y <- trainY_type2
    gamma_0 <- type1_gamma_0
    gamma_1 <- type1_gamma_1
    prior <- type1_prior

orig_X <- X
orig_Y <- Y

recom_list <- c()
recom_instant_prob <- c()
for (pep_idx in 1:no_recom) {
    print (pep_idx)
    print (which_type)
    print ("add_ins")
    print (add_ins)
    print ("mean theta")

    alpha <- Dirichlet_Parameter(X, Y, AAclass, gamma_0, gamma_1)

    # generate random peptide lib
    pep_lib <- gen_peptide_lib(Nlib, nAA, nL+nR, maxL, maxR, minL, minR, S.Pos,colnames(X))
    print ("random peptide generated")

    theta <- getTheta(alpha = alpha, classlist = AAclass)
    prob <- NB_predict(pep_lib, theta, S.Pos, maxL, maxR, prior)

    pep_list <- select_new_recom(X, pep_lib, prob)
    print ("selected new recom")
    new_pep <- pep_list$peptide
    new_pep_prob <- pep_list$prob
    print (new_pep_prob)
    recom_list <- rbind(recom_list, new_pep)
    recom_instant_prob <- c(recom_instant_prob, new_pep_prob)
    # add to training data
    for (add in 1:add_ins) {
        X <- rbind(X, new_pep)
        Y <- c(Y, 0)
    }
}

# validate prob of recom_list
#Parallel Computation Packages
require(multicore)
require(snow)
require(doParallel)
#=================================================================================
#Parallel Computation Setting
# hosts <- rep('whale',40)
# hosts <- c(rep("whale",52),rep("ahab",24),rep("flask",24),rep("starbuck",24))
# hosts <- c(rep("whale",30),rep("ahab",30),rep("flask",30),rep("starbuck",30))
# hosts <- c(rep("whale",10),rep("ahab",10),rep("flask",10),rep("starbuck",10))
# hosts <- c(rep("ishmael",30),rep("stubb",30),rep("daggoo",30),rep("tashtego",30))
hosts <- c(rep("ishmael",10),rep("stubb",10),rep("daggoo",10),rep("tashtego",10))
cl <- makeCluster(hosts, type = "SOCK")
registerDoParallel(cl)

iteration <- 500
alpha <- Dirichlet_Parameter(orig_X, orig_Y, AAclass, gamma_0, gamma_1)

colnames(recom_list) <- colnames(orig_X)

predict_mat <- foreach (n = 1:iteration, .init = c(), .combine = rbind, .inorder=F, .errorhandling='pass') %dopar% {
    theta <- getTheta_MC(alpha = alpha, classlist = AAclass)
    NB_predict(recom_list, theta, S.Pos, maxL, maxR, prior)
}
recom_prob <- colMeans(predict_mat)

stopCluster(cl)
