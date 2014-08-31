rm(list=ls())

root_path <- "/fs/home/jw865/peptide-catalysis"
load(paste(root_path, "/data/all_data_for_training.RData", sep=''))

# parameters that can change
no_recom <- 101

sfp_gamma_0 <- 100
sfp_gamma_1 <- 0.5
sfp_prior <- 1e-4

AcpS_gamma_0 <- 100
AcpS_gamma_1 <- 0.5
AcpS_prior <- 1e-4

PfAcpH_gamma_0 <- 10
PfAcpH_gamma_1 <- 0.05
PfAcpH_prior <- 0.5

add_ins <- 5
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

nAA <- length(unique(as.numeric(AAclass)))
X_unlabel <- trainX_PfAcpH
Y_unlabel <- trainY_PfAcpH
unlabel_gamma_0 <- PfAcpH_gamma_0
unlabel_gamma_1 <- PfAcpH_gamma_1
unlabel_prior <- PfAcpH_prior
# find type 1 hit
    which_type <- 'type 1 hit'
    X_label_prefer <- trainX_sfp
    Y_label_prefer <- trainY_sfp
    X_label_unprefer <- trainX_AcpS
    Y_label_unprefer <- trainY_AcpS
    label_prefer_gamma_0 <- sfp_gamma_0
    label_prefer_gamma_1 <- sfp_gamma_1
    label_prefer_prior <- sfp_prior
    label_unprefer_gamma_0 <- AcpS_gamma_0
    label_unprefer_gamma_1 <- AcpS_gamma_1
    label_unprefer_prior <- AcpS_prior
# find type 2 hit
#     which_type <- 'type 2 hit'
#     X_label_prefer <- trainX_AcpS
#     Y_label_prefer <- trainY_AcpS
#     X_label_unprefer <- trainX_sfp
#     Y_label_unprefer <- trainY_sfp
#     label_prefer_gamma_0 <- AcpS_gamma_0
#     label_prefer_gamma_1 <- AcpS_gamma_1
#     label_prefer_prior <- AcpS_prior
#     label_unprefer_gamma_0 <- sfp_gamma_0
#     label_unprefer_gamma_1 <- sfp_gamma_1
#     label_unprefer_prior <- sfp_prior

orig_X_label_prefer <- X_label_prefer
orig_Y_label_prefer <- Y_label_prefer
orig_X_label_unprefer <- X_label_unprefer
orig_Y_label_unprefer <- Y_label_unprefer
orig_X_unlabel <- X_unlabel
orig_Y_unlabel <- Y_unlabel

recom_list <- c()
recom_instant_prob <- c()
for (pep_idx in 1:no_recom) {
    print (pep_idx)
    print (which_type)
    print ("add_ins")
    print (add_ins)
    print ("mean theta")

    alpha_label_prefer <- Dirichlet_Parameter(X_label_prefer, Y_label_prefer, AAclass, label_prefer_gamma_0, label_prefer_gamma_1)
    alpha_label_unprefer <- Dirichlet_Parameter(X_label_unprefer, Y_label_unprefer, AAclass, label_unprefer_gamma_0, label_unprefer_gamma_1)
    alpha_unlabel <- Dirichlet_Parameter(X_unlabel, Y_unlabel, AAclass, unlabel_gamma_0, unlabel_gamma_1)

    # generate random peptide lib
    pep_lib <- gen_peptide_lib(Nlib, nAA, nL+nR, maxL, maxR, minL, minR, S.Pos,colnames(X_label_prefer))
    print ("random peptide generated")

    label_prefer_theta <- getTheta(alpha = alpha_label_prefer, classlist = AAclass)
    label_unprefer_theta <- getTheta(alpha = alpha_label_unprefer, classlist = AAclass)
    unlabel_theta <- getTheta(alpha = alpha_unlabel, classlist = AAclass)
    prob <- new_NB_predict(pep_lib, label_prefer_theta, label_unprefer_theta, unlabel_theta, S.Pos, maxL, maxR, label_prefer_prior, label_unprefer_prior, unlabel_prior)

    pep_list <- select_new_recom(rbind(X_label_prefer, X_label_unprefer, X_unlabel), pep_lib, prob)
    print ("selected new recom")
    new_pep <- pep_list$peptide
    new_pep_prob <- pep_list$prob
    print (new_pep_prob)
    recom_list <- rbind(recom_list, new_pep)
    recom_instant_prob <- c(recom_instant_prob, new_pep_prob)
    # add to training data
    for (add in 1:add_ins) {
        X_label_prefer <- rbind(X_label_prefer, new_pep)
        Y_label_prefer <- c(Y_label_prefer, 0)
        X_label_unprefer <- rbind(X_label_unprefer, new_pep)
        Y_label_unprefer <- c(Y_label_unprefer, 0)
        X_unlabel <- rbind(X_unlabel, new_pep)
        Y_unlabel <- c(Y_unlabel, 0)
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
hosts <- c(rep("whale",10),rep("ahab",10),rep("flask",10),rep("starbuck",10))
# hosts <- c(rep("ishmael",30),rep("stubb",30),rep("daggoo",30),rep("tashtego",30))
# hosts <- c(rep("ishmael",10),rep("stubb",10),rep("daggoo",10),rep("tashtego",10))
cl <- makeCluster(hosts, type = "SOCK")
registerDoParallel(cl)

iteration <- 500
alpha_label_prefer <- Dirichlet_Parameter(orig_X_label_prefer, orig_Y_label_prefer, AAclass, label_prefer_gamma_0, label_prefer_gamma_1)
alpha_label_unprefer <- Dirichlet_Parameter(orig_X_label_unprefer, orig_Y_label_unprefer, AAclass, label_unprefer_gamma_0, label_unprefer_gamma_1)
alpha_unlabel <- Dirichlet_Parameter(orig_X_unlabel, orig_Y_unlabel, AAclass, unlabel_gamma_0, unlabel_gamma_1)

colnames(recom_list) <- colnames(orig_X_label_prefer)

predict_mat <- foreach (n = 1:iteration, .init = c(), .combine = rbind, .inorder=F, .errorhandling='pass') %dopar% {
    label_prefer_theta <- getTheta_MC(alpha = alpha_label_prefer, classlist = AAclass)
    label_unprefer_theta <- getTheta_MC(alpha = alpha_label_unprefer, classlist = AAclass)
    unlabel_theta <- getTheta_MC(alpha = alpha_unlabel, classlist = AAclass)
    new_NB_predict(recom_list, label_prefer_theta, label_unprefer_theta, unlabel_theta, S.Pos, maxL, maxR, label_prefer_prior, label_unprefer_prior, unlabel_prior)
}
recom_prob <- colMeans(predict_mat)

stopCluster(cl)
