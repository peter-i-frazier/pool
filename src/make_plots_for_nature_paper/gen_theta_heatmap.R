rm(list = ls())
root <- "/fs/home/jw865/remote_deployment/ucsd_reversible_labeling/src/"
source(paste0(root, "core/dependency.R"))
ResolveDependency(root)

label_name <- 'sfp'
notlabel_name <- 'AcpS'
num.mc = 10000
cnames = c('L19', 'L18', 'L17', 'L16', 'L15', 'L14', 'L13', 'L12', 'L11', 'L10', 'L9', 'L8', 'L7', 'L6', 'L5', 'L4', 'L3', 'L2', 'L1',
'R1', 'R2', 'R3', 'R4', 'R5', 'R6', 'R7', 'R8', 'R9', 'R10', 'R11', 'R12', 'R13', 'R14', 'R15', 'R16', 'R17', 'R18', 'R19')

#prior.df <- data.frame(sfp = c(10, 0.1, 1e-4), AcpS = c(10, 1, 1e-4), PfAcpH = c(10, 10, 0.5))
prior.df <- data.frame(sfp = c(1, 1, 1e-4), AcpS = c(1, 1, 1e-4), PfAcpH = c(10, 10, 0.5))
label.alpha.0 <- SetPriorReducedAA(prior.df[1, label_name], NUM_CLASS)
label.alpha.1 <- SetPriorReducedAA(prior.df[2, label_name], NUM_CLASS)
label.p1 <- prior.df[3, label_name]

notlabel.alpha.0 <- SetPriorReducedAA(prior.df[1, notlabel_name], NUM_CLASS)
notlabel.alpha.1 <- SetPriorReducedAA(prior.df[2, notlabel_name], NUM_CLASS)
notlabel.p1 <- prior.df[3, notlabel_name]

true_dataset <- GetDataReducedAA("SELECT * FROM binary_labeling_activity")
true_label_idx <- true_dataset$data[, label_name] != -1
true_notlabel_idx <- true_dataset$data[, notlabel_name] != -1
true_PfAcpH_idx <- true_dataset$data[, 'PfAcpH'] != -1
true_label_X = true_dataset$feature[true_label_idx, ]
true_label_Y = true_dataset$data[true_label_idx, label_name]
true_notlabel_X = true_dataset$feature[true_notlabel_idx, ]
true_notlabel_Y = true_dataset$data[true_notlabel_idx, notlabel_name]

label.trained.params <- BayesianNaiveBayes(true_label_X, true_label_Y, label.alpha.1, label.alpha.0, label.p1)
label.thetas.1 <- SampleThetas(label.trained.params$post.alpha.1, num.mc)
label.thetas.0 <- SampleThetas(label.trained.params$post.alpha.0, num.mc)

notlabel.trained.params <- BayesianNaiveBayes(true_notlabel_X, true_notlabel_Y, notlabel.alpha.1, notlabel.alpha.0, notlabel.p1)
notlabel.thetas.1 <- SampleThetas(notlabel.trained.params$post.alpha.1, num.mc)
notlabel.thetas.0 <- SampleThetas(notlabel.trained.params$post.alpha.0, num.mc)

result_mean = c()
result_sd = c()
label_theta_1 = c()
label_theta_0 = c()
label_theta_1_sd = c()
label_theta_0_sd = c()
not_label_theta_0 = c()
not_label_theta_1 = c()
not_label_theta_0_sd = c()
not_label_theta_1_sd = c()
label_ratio_mean = c()
label_ratio_sd = c()
notlabel_ratio_mean = c()
notlabel_ratio_sd = c()
label_alpha_1 = c()
label_alpha_0 = c()
notlabel_alpha_1 = c()
notlabel_alpha_0 = c()
for (k in 1:38) {
  if (k < 20) {
    j = 20-k
  } else {
    j = k
  }
  label_alpha_1 = cbind(label_alpha_1, label.trained.params$post.alpha.1[[j]])
  label_alpha_0 = cbind(label_alpha_0, label.trained.params$post.alpha.0[[j]])
  notlabel_alpha_1 = cbind(notlabel_alpha_1, notlabel.trained.params$post.alpha.1[[j]])
  notlabel_alpha_0 = cbind(notlabel_alpha_0, notlabel.trained.params$post.alpha.0[[j]])

  label_theta_1 = cbind(label_theta_1, apply(label.thetas.1[[j]], 2, mean))
  label_theta_0 = cbind(label_theta_0, apply(label.thetas.0[[j]], 2, mean))
  label_theta_1_sd = cbind(label_theta_1_sd, apply(label.thetas.1[[j]], 2, sd))
  label_theta_0_sd = cbind(label_theta_0_sd, apply(label.thetas.0[[j]], 2, sd))
  not_label_theta_1 = cbind(not_label_theta_1, apply(notlabel.thetas.1[[j]], 2, mean))
  not_label_theta_0 = cbind(not_label_theta_0, apply(notlabel.thetas.0[[j]], 2, mean))
  not_label_theta_1_sd = cbind(not_label_theta_1_sd, apply(notlabel.thetas.1[[j]], 2, sd))
  not_label_theta_0_sd = cbind(not_label_theta_0_sd, apply(notlabel.thetas.0[[j]], 2, sd))
  ratio = log(label.thetas.1[[j]]) - log(label.thetas.0[[j]]) - log(notlabel.thetas.1[[j]]) + log(notlabel.thetas.0[[j]])
  ratio_mean = apply(ratio, 2, mean)
  #ratio_mean = ratio_mean - mean(ratio_mean)
  ratio_sd = apply(ratio, 2, sd)
  result_mean = cbind(result_mean, ratio_mean)
  result_sd = cbind(result_sd, ratio_sd)
  label_ratio = log(label.thetas.1[[j]]) - log(label.thetas.0[[j]])
  label_ratio_mean = cbind(label_ratio_mean, apply(label_ratio, 2, mean))
  label_ratio_sd = cbind(label_ratio_sd, apply(label_ratio, 2, sd))
  notlabel_ratio = log(notlabel.thetas.1[[j]]) - log(notlabel.thetas.0[[j]])
  notlabel_ratio_mean = cbind(notlabel_ratio_mean, apply(notlabel_ratio, 2, mean))
  notlabel_ratio_sd = cbind(notlabel_ratio_sd, apply(notlabel_ratio, 2, sd))
}
colnames(result_mean) <- cnames
colnames(result_sd) <- cnames
colnames(label_alpha_1) = cnames
colnames(label_alpha_0) = cnames
colnames(notlabel_alpha_1) = cnames
colnames(notlabel_alpha_0) = cnames

write.csv(label_alpha_1, "label_alpha_1.csv")
write.csv(label_alpha_0, "label_alpha_0.csv")
write.csv(notlabel_alpha_1, "notlabel_alpha_1.csv")
write.csv(notlabel_alpha_0, "notlabel_alpha_0.csv")

#write.csv(result_mean, "heatmap_mean.csv")
#write.csv(result_sd, "heatmap_sd.csv")
#table <- result_mean / result_sd
#colnames(table) <- cnames
#write.csv(table, "heatmap.csv")
#write.csv(label_theta_1, "label_theta_1.csv")
#write.csv(label_theta_0, "label_theta_0.csv")
#write.csv(label_theta_1_sd, "label_theta_1_sd.csv")
#write.csv(label_theta_0_sd, "label_theta_0_sd.csv")
#write.csv(not_label_theta_1, "not_label_theta_1.csv")
#write.csv(not_label_theta_0, "not_label_theta_0.csv")
#write.csv(not_label_theta_0_sd, "not_label_theta_0_sd.csv")
#write.csv(not_label_theta_1_sd, "not_label_theta_1_sd.csv")
#write.csv(notlabel_ratio_mean, "notlabel_logratio_mean.csv")
#write.csv(notlabel_ratio_sd, "notlabel_logratio_sd.csv")
#write.csv(label_ratio_mean, "label_logratio_mean.csv")
#write.csv(label_ratio_sd, "label_logratio_sd.csv")
