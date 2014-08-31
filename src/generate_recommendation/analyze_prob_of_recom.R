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
# cl <- makeCluster(8)
registerDoParallel(cl)

colnames(recom_list) <- colnames(trainX_sfp)

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
#    which_type <- 'type 2 hit'
#    X_label_prefer <- trainX_AcpS
#    Y_label_prefer <- trainY_AcpS
#    X_label_unprefer <- trainX_sfp
#    Y_label_unprefer <- trainY_sfp
#    label_prefer_gamma_0 <- AcpS_gamma_0
#    label_prefer_gamma_1 <- AcpS_gamma_1
#    label_prefer_prior <- AcpS_prior
#    label_unprefer_gamma_0 <- sfp_gamma_0
#    label_unprefer_gamma_1 <- sfp_gamma_1
#    label_unprefer_prior <- sfp_prior

alpha_label_prefer <- Dirichlet_Parameter(X_label_prefer, Y_label_prefer, AAclass, label_prefer_gamma_0, label_prefer_gamma_1)
alpha_label_unprefer <- Dirichlet_Parameter(X_label_unprefer, Y_label_unprefer, AAclass, label_unprefer_gamma_0, label_unprefer_gamma_1)
alpha_unlabel <- Dirichlet_Parameter(X_unlabel, Y_unlabel, AAclass, unlabel_gamma_0, unlabel_gamma_1)


# label_prefer
label_prefer_predict_mat <- foreach (n = 1:iteration, .init = c(), .combine = rbind, .inorder=F, .errorhandling='pass') %dopar% {
    theta <- getTheta_MC(alpha = alpha_label_prefer, classlist = AAclass)
    NB_predict(recom_list, theta, S.Pos, maxL, maxR, label_prefer_prior)
}
label_prefer_predProb <- colMeans(label_prefer_predict_mat)
# label_unprefer
label_unprefer_predict_mat <- foreach (n = 1:iteration, .init = c(), .combine = rbind, .inorder=F, .errorhandling='pass') %dopar% {
    theta <- getTheta_MC(alpha = alpha_label_unprefer, classlist = AAclass)
    NB_predict(recom_list, theta, S.Pos, maxL, maxR, label_unprefer_prior)
}
label_unprefer_predProb <- colMeans(label_unprefer_predict_mat)
# unlabel
unlabel_predict_mat <- foreach (n = 1:iteration, .init = c(), .combine = rbind, .inorder=F, .errorhandling='pass') %dopar% {
    theta <- getTheta_MC(alpha = alpha_unlabel, classlist = AAclass)
    NB_predict(recom_list, theta, S.Pos, maxL, maxR, unlabel_prior)
}
unlabel_predProb <- colMeans(unlabel_predict_mat)

recom_prob <- label_prefer_predProb * (1-label_unprefer_predProb) * unlabel_predProb
