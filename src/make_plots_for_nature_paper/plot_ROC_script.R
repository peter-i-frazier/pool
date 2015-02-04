root_path = "/fs/home/jw865/peptide-catalysis"
datafile <- paste(root_path,"/data/all_TS_data.csv",sep='')

gamma_0_sfp <- 100
gamma_1_sfp <- 0.01
prior_sfp <- 1e-4

gamma_0_AcpS <- 100
gamma_1_AcpS <- 0.1
prior_AcpS <- 1e-4

gamma_0_PfAcpH <- 1000
gamma_1_PfAcpH <- 10
prior_PfAcpH <- 1e-4

# gamma_0_type1 <- 50
# gamma_1_type1 <- 0.5
# prior_type1 <- 1e-4
# 
# gamma_0_type2 <- 1000
# gamma_1_type2 <- 0.5
# prior_type2 <- 1e-4

nL <- 19
nR <- 19
S.Pos <- nL
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

X_PfAcpH <- train_data[train_data[,'PfAcpH'] != -1,1:(nL+nR)]
Y_PfAcpH <- train_data[train_data[,'PfAcpH'] != -1, 'PfAcpH']
X_sfp <- train_data[train_data[,'sfp'] != -1, 1:(nL+nR)]
Y_sfp <- train_data[train_data[,'sfp'] != -1, 'sfp']
X_AcpS <- train_data[train_data[,'AcpS'] != -1, 1:(nL+nR)]
Y_AcpS <- train_data[train_data[,'AcpS'] != -1, 'AcpS']
X_type1 <- train_data[train_data[,'type1'] != -1, 1:(nL+nR)]
Y_type1 <- train_data[train_data[,'type1'] != -1, 'type1']
X_type2 <- train_data[train_data[,'type2'] != -1, 1:(nL+nR)]
Y_type2 <- train_data[train_data[,'type2'] != -1, 'type2']


## type1 cv
#type1_cv <- loocv(X_type1, Y_type1, AAclass, S.Pos, nL, nR, gamma_0_type1, gamma_1_type1, prior_type1, itr) 
#xy <- type1_cv$axis
#pdf("type1_hit_ROC.pdf")
#plot(x = xy$x, y = xy$y, type = 'l', xlim = c(0,1), ylim = c(0,1), xlab = "false positive rate", ylab = "true positive rate")
#title(main = "type1 hit ROC")
#dev.off()

# # type2 cv
# type2_cv <- loocv(X_type2, Y_type2, AAclass, S.Pos, nL, nR, gamma_0_type2, gamma_1_type2, prior_type2, itr) 
# xy <- type2_cv$axis
# pdf("type2_hit_ROC.pdf")
# plot(x = xy$x, y = xy$y, type = 'l', xlim = c(0,1), ylim = c(0,1), xlab = "false positive rate", ylab = "true positive rate")
# title(main = "type2 hit ROC")
# dev.off()

# new method type 1 ROC
# prob_sfp <- Naive_Bayes(X_sfp, Y_sfp, X_type1, AAclass, S.Pos, nL, nR, gamma_0_sfp, gamma_1_sfp, prior_sfp, itr)
# prob_AcpS <- Naive_Bayes(X_AcpS, Y_AcpS, X_type1, AAclass, S.Pos, nL, nR, gamma_0_AcpS, gamma_1_AcpS, prior_AcpS, itr)
# prob_PfAcpH<- Naive_Bayes(X_PfAcpH, Y_PfAcpH, X_type1, AAclass, S.Pos, nL, nR, gamma_0_PfAcpH, gamma_1_PfAcpH, prior_PfAcpH, itr)
# prob <- prob_sfp * (1. - prob_AcpS) * prob_PfAcpH
# xy <- ROC_xy_output(prob, Y_type1)
# pdf("type1_hit_new_method_ROC.pdf")
# plot(x = xy$x, y = xy$y, type = 'l', xlim = c(0,1), ylim = c(0,1), xlab = "false positive rate", ylab = "true positive rate")
# title(main = "type1 hit new method ROC")
# dev.off()
# table <- cbind(xy$x, xy$y)
# colnames(table) <- c("type1_new_x", "type1_new_y")
# write.csv(table, "data_for_plot_type1.csv")
# 
# new method type 2 ROC
prob_sfp <- Naive_Bayes(X_sfp, Y_sfp, X_type2, AAclass, S.Pos, nL, nR, gamma_0_sfp, gamma_1_sfp, prior_sfp, itr)
prob_AcpS <- Naive_Bayes(X_AcpS, Y_AcpS, X_type2, AAclass, S.Pos, nL, nR, gamma_0_AcpS, gamma_1_AcpS, prior_AcpS, itr)
prob_PfAcpH<- Naive_Bayes(X_PfAcpH, Y_PfAcpH, X_type2, AAclass, S.Pos, nL, nR, gamma_0_PfAcpH, gamma_1_PfAcpH, prior_PfAcpH, itr)
prob <- prob_AcpS * (1. - prob_sfp) * prob_PfAcpH
xy <- ROC_xy_output(prob, Y_type2)
pdf("type2_hit_new_method_ROC.pdf")
plot(x = xy$x, y = xy$y, type = 'l', xlim = c(0,1), ylim = c(0,1), xlab = "false positive rate", ylab = "true positive rate")
title(main = "type2 hit new method ROC")
dev.off()
table <- cbind(xy$x, xy$y)
colnames(table) <- c("type2_new_x", "type2_new_y")
write.csv(table, "data_for_plot_type2.csv")

# new method sfp
# print ("plot sfp")
# num_fold <- 260
# cv <- fold_cv(X_sfp, Y_sfp, AAclass, S.Pos, nL, nR, gamma_0_sfp, gamma_1_sfp, prior_sfp, itr, num_fold)
# xy <- cv$xy
# pdf("sfp_new_method_ROC.pdf")
# plot(x = xy$x, y = xy$y, type = 'l', xlim = c(0,1), ylim = c(0,1), xlab = "false positive rate", ylab = "true positive rate")
# title(main = "sfp new method ROC")
# dev.off()
# table <- cbind(xy$x, xy$y)
# colnames(table) <- c("sfp_new_x", "sfp_new_y")
# write.csv(table, "data_for_plot_sfp.csv")
# 
# new method AcpS
# print ("plot AcpS")
# num_fold <- 260
# cv <- fold_cv(X_AcpS, Y_AcpS, AAclass, S.Pos, nL, nR, gamma_0_AcpS, gamma_1_AcpS, prior_AcpS, itr, num_fold)
# xy <- cv$xy
# pdf("AcpS_new_method_ROC.pdf")
# plot(x = xy$x, y = xy$y, type = 'l', xlim = c(0,1), ylim = c(0,1), xlab = "false positive rate", ylab = "true positive rate")
# title(main = "AcpS new method ROC")
# dev.off()
# table <- cbind(xy$x, xy$y)
# colnames(table) <- c("AcpS_new_x", "AcpS_new_y")
# write.csv(table, "data_for_plot_AcpS.csv")

# # new method PfAcpH
# print ("plot PfAcpH")
# num_fold <- 500
# cv <- fold_cv(X_PfAcpH, Y_PfAcpH, AAclass, S.Pos, nL, nR, gamma_0_PfAcpH, gamma_1_PfAcpH, prior_PfAcpH, itr, num_fold)
# xy <- cv$xy
# pdf("PfAcpH_new_method_ROC.pdf")
# plot(x = xy$x, y = xy$y, type = 'l', xlim = c(0,1), ylim = c(0,1), xlab = "false positive rate", ylab = "true positive rate")
# title(main = "PfAcpH new method ROC")
# dev.off()
# table <- cbind(xy$x, xy$y)
# colnames(table) <- c("PfAcpH_new_x", "PfAcpH_new_y")
# write.csv(table, "data_for_plot_PfAcpH.csv")
