rm(list=ls())
source("normalization_util.R")
source("simulate_data.R")

toy_prob <- function(data, j_class, treatment) {
# Construct matrix D and A for quadratic program with linear constraint, which has form
# min x^T D x, s.t A^T x >= 0
    y_index <- matrix(rep(-1, 2*dim(data)[1]), nrow=dim(data)[1], ncol=2)
    colnames(y_index) <- c('i', 'j')
    unique_seqs <- c()
    for (n in 1:dim(data)[1]) {
        if (return_idx(data[n, 'seq'], unique_seqs) == 0) {
            unique_seqs <- c(unique_seqs, data[n, 'seq'])
        }
        y_index[n, 'i'] <- return_idx(data[n, 'seq'], unique_seqs)
        # index j is not really TS, it's the index of TS in j_class
        y_index[n, 'j'] <- which(j_class == data[n, 'TS'])
    }
    num_i <- length(unique_seqs)
    num_j <- length(j_class)
    len_x <- 2 * num_j + num_i
    # construct D matrix
    D_matrix <- matrix(rep(0, len_x * len_x), nrow=len_x, ncol=len_x)
    for (n in 1:dim(data)[1]) {
        idx_alpha <- y_index[n, 'j']
        idx_beta <- num_j + y_index[n, 'j']
        idx_theta <- 2 * num_j + y_index[n, 'i']
        D_matrix[idx_alpha, idx_alpha] <- data[n, treatment]^2 + D_matrix[idx_alpha, idx_alpha]
        D_matrix[idx_beta, idx_beta] <- 1.0 + D_matrix[idx_beta, idx_beta]
        D_matrix[idx_theta, idx_theta] <- 1.0 + D_matrix[idx_theta, idx_theta]
        D_matrix[idx_alpha, idx_beta] <- -data[n, treatment] + D_matrix[idx_alpha, idx_beta]
        D_matrix[idx_beta, idx_alpha] <- -data[n, treatment] + D_matrix[idx_beta, idx_alpha]
        D_matrix[idx_alpha, idx_theta] <- -data[n, treatment] + D_matrix[idx_alpha, idx_theta]
        D_matrix[idx_theta, idx_alpha] <- -data[n, treatment] + D_matrix[idx_theta, idx_alpha]
        D_matrix[idx_beta, idx_theta] <- 1.0 + D_matrix[idx_beta, idx_theta]
        D_matrix[idx_theta, idx_beta] <- 1.0 + D_matrix[idx_theta, idx_beta]
    }
    # construct A matrix: A_matrix %*% x >= b_0, where first meq are equality constraints
    # equality part
    # # sum of theta = constant
    A_sum <- matrix(c(rep(0, 2 * num_j), rep(1, num_i)), nrow=1)
    meq <- 1
    A_alpha_beta <- matrix(rep(0, 2 * num_j * len_x), nrow=2*num_j)
    for (i in 1:(2*num_j)) {
        A_alpha_beta[i, i] <- 1
    }
    #A_matrix <- rbind(A_sum, A_alpha_beta)
    A_matrix <- t(A_sum)
    # b_0
    b_0 <- c(num_i, rep(0, 2*num_j))
    # dvec
    dvec <- matrix(rep(0, 2*num_j + num_i), nrow=1)
    return (list(D=D_matrix, A=t(A_matrix), b_0=b_0, unique_seqs=unique_seqs, meq=meq, dvec=t(dvec)))
}

# quadratic approach
require(kernlab)
root_path = "../.."
org.data <- read.csv(paste(root_path, "/data/all_TS_raw_reading.csv", sep=''), stringsAsFactors=F)
# unique_seqs <- toy_prob(org.data, 1:5, 'sfp_1')$unique_seqs
# dup_data <- c()
# for (seq in unique_seqs) {
#     sub_data <- org.data[org.data[,'seq'] == seq,]
#     if (dim(sub_data)[1] > 1) {
#         dup_data <- rbind(dup_data, sub_data)
#     }
# }
# write.csv(dup_data, 'dup_data.csv', row.names=F)
dup_data <- read.csv('dup_data.csv', stringsAsFactors=F)

treatment <- "sfp_1"
j_class <- c(1,2,3,4,5)
#j_class <- c(1,2)
eps <- 1e-3
bound <- 1e7
true_sigma <- 1:5
true_mu <- c(0, 10 * 2:5)
num_data <- 30
sum_theta <- simulate_data(true_mu, true_sigma, num_data, j_class)

dataset <- dup_data[!is.na(dup_data[, treatment]),]
dataset[dataset[,'TS']==3,treatment] <- dataset[dataset[,'TS']==3,treatment] + 10000
#dataset <- read.csv('simulate_data.csv', stringsAsFactors=F)
matrix_list <- toy_prob(dataset, j_class, treatment)
D <- matrix_list$D
A <- matrix_list$A
unique_seqs <- matrix_list$unique_seqs
ll <- dim(D)[1]
num_j <- length(j_class)
num_i <- ll - 2*num_j
small_val <- 0.01
# ### find no_label idx
# nolabel_table <- read.csv(paste(root_path, "/data/nolabel_table.csv", sep=''), stringsAsFactors=F)
# nolabel_table <- nolabel_table[nolabel_table[,'type'] == treatment,]
# nolabel_list <- rep(0, num_i)
# for (n in 1:dim(nolabel_table)[1]) {
#     seq <- org.data[org.data[,'TS'] == nolabel_table[n, 'TS'] & org.data[,'spot'] == paste(nolabel_table[n, 'row'], as.character(nolabel_table[n, 'col']), sep=''), 'seq']
#     nolabel_list[return_idx(seq, unique_seqs)] <- 1
# }
# lower_bound <- c(1, rep(0, num_j-1), 0, rep(-bound, num_j+num_i-1))
# upper_bound <- c(1+eps, rep(bound, num_j-1), eps, rep(bound, num_j+num_i-1))
lower_bound <- c(rep(0,2*num_j), rep(-bound, num_i))
upper_bound <- rep(bound, 2*num_j+num_i)
A <- matrix(0, nrow=2, ncol=ll)
A[1,1] <- 1
A[2,6] <- 1
bvec <- c(1, 0)
ss <- LowRankQP(Vmat=D, dvec=rep(0,ll), Amat=A, bvec=bvec, uvec=upper_bound, method="LU", verbose=F)
print (1/ss$alpha[1:5])
print (ss$alpha[6:10]/ss$alpha[1:5])
print (t(ss$alpha) %*% D %*% ss$alpha)
#for (i in 1:num_i) {
#    if (nolabel_list[i] == 1) {
#        upper_bound[2*num_j + i] <- small_val
#    }
#}
# ##### TEST
# A <- matrix(rep(0,ll*3), nrow=3)
# A[1,1] <- 1
# A[2,2] <- 1
# A[3,3] <- 1
# #####
# perturb_mat <- eps * diag(ll)
# for (i in 1:(2*num_j)) { perturb_mat[i,i] <- 0 }
# H <- D+perturb_mat
# sol <- ipop(c=rep(0,ll), H=H, A=A, b=c(1,0,0), l=lower_bound, u=upper_bound, r=c(1e-6,bound,1e-6), margin=0.05, verb=1)
# theta_vec <- primal(sol)[(2 * num_j + 1):ll]
# sigma_vec <- 1 / primal(sol)[1:num_j]
# mu_vec <- sigma_vec * primal(sol)[(num_j + 1):(2 * num_j)]
# print ('mu_vec')
# print (mu_vec)
# print ('sigma_vec')
# print (sigma_vec)
# # true vec
# seed <- read.csv('simulate_seed.csv', stringsAsFactors=F)
# true_theta <- seed[, 'sfp_1'][1:num_i]
# true_vec <- c(1/true_sigma, true_mu/true_sigma, true_theta)
# # print obj value
# print ('optim vs true on D')
# print (0.5 * t(primal(sol)) %*% D %*% primal(sol))
# print (0.5 * t(true_vec) %*% D %*% true_vec)
# print ('optim vs true on H')
# print (0.5 * t(primal(sol)) %*% H %*% primal(sol))
# print (0.5 * t(true_vec) %*% H %*% true_vec)
# 
# 
# # construct norm result table
# norm.data <- cbind(dataset[,c('TS', 'spot', 'seq')], matrix(rep(NA, 6 * dim(dataset)[1]), ncol=6))
# colnames(norm.data) <- c('TS', 'spot', 'seq', 'unique_seq_idx', 'raw', 'theta', 'eps', 'mu', 'sigma')
# for (n in 1:dim(dataset)[1]) {
#     idx <- return_idx(norm.data[n, 'seq'], unique_seqs)
#     norm.data[n, 'raw'] <- dataset[n, treatment]
#     norm.data[n, 'unique_seq_idx'] <- idx
#     norm.data[n, 'theta'] <- theta_vec[idx]
#     norm.data[n, 'sigma'] <- sigma_vec[which(j_class == norm.data[n, 'TS'])]
#     norm.data[n, 'mu'] <- mu_vec[which(j_class == norm.data[n, 'TS'])]
#     norm.data[n, 'eps'] <- (norm.data[n, 'raw'] - norm.data[n, 'mu']) / norm.data[n, 'sigma'] - norm.data[n, 'theta']
# }
# write.csv(norm.data, "toy_norm_result_table.csv", row.names=F)
