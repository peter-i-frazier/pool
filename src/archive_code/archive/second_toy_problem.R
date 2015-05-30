rm(list=ls())
source("normalization_util.R")
source("simulate_data.R")

toy_prob <- function(data, j_class, j_zero, treatment) {
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
        if (data[n, 'TS'] == j_zero) {
            y_index[n, 'j'] <- 0
        } else {
            y_index[n, 'j'] <- which(j_class == data[n, 'TS'])
        }
    }
    num_i <- length(unique_seqs)
    num_j <- length(j_class)
    len_x <- 2 * num_j + num_i
    # construct D matrix and dvec
    D_matrix <- matrix(rep(0, len_x * len_x), nrow=len_x, ncol=len_x)
    dvec <- rep(0, len_x)
    sum_y02 <- 0
    for (n in 1:dim(data)[1]) {
        idx_theta <- 2 * num_j + y_index[n, 'i']
        D_matrix[idx_theta, idx_theta] <- 1.0 + D_matrix[idx_theta, idx_theta]
        if (y_index[n, 'j'] > 0) {
            idx_alpha <- y_index[n, 'j']
            idx_beta <- num_j + y_index[n, 'j']
            D_matrix[idx_alpha, idx_alpha] <- data[n, treatment]^2 + D_matrix[idx_alpha, idx_alpha]
            D_matrix[idx_beta, idx_beta] <- 1.0 + D_matrix[idx_beta, idx_beta]
            D_matrix[idx_alpha, idx_beta] <- -data[n, treatment] + D_matrix[idx_alpha, idx_beta]
            D_matrix[idx_beta, idx_alpha] <- -data[n, treatment] + D_matrix[idx_beta, idx_alpha]
            D_matrix[idx_alpha, idx_theta] <- -data[n, treatment] + D_matrix[idx_alpha, idx_theta]
            D_matrix[idx_theta, idx_alpha] <- -data[n, treatment] + D_matrix[idx_theta, idx_alpha]
            D_matrix[idx_beta, idx_theta] <- 1.0 + D_matrix[idx_beta, idx_theta]
            D_matrix[idx_theta, idx_beta] <- 1.0 + D_matrix[idx_theta, idx_beta]
        } else {
            dvec[idx_theta] <- dvec[idx_theta] + 2 * data[n, treatment]
            sum_y02 <- sum_y02 + data[n, treatment]^2
        }
    }
    tA <- cbind(diag(num_j), matrix(0, nrow=(num_j), ncol=(num_j+num_i)))
    b_0 <- rep(0, num_j)
    return (list(D=2*D_matrix, dvec=dvec, A=t(tA), b_0=b_0, sum_y02=sum_y02, unique_seqs=unique_seqs))
}

# quadratic approach
require(quadprog)
require(modeest)
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
j_class <- c(2,3,4,5)
eps <- 1e-3
bound <- 1e7
true_sigma <- c(5, 3, 2, 10, 9)
true_mu <- c(15, 10, 1, 9, 20)
true_mu_vec <- true_mu[2:5] - true_sigma[2:5] / true_sigma[1] * true_mu[1]
true_sigma_vec <- true_sigma[2:5] / true_sigma[1]
num_data <- 10
true_theta <- simulate_data(true_mu, true_sigma, num_data, 1:5)

dataset <- org.data[!is.na(org.data[, treatment]),]
for (ts in 1:5) {
    mode <- mlv(dataset[dataset[,'TS']==ts,treatment], method='lientz', bw=0.2)$M
    std <- sd(dataset[dataset[,'TS']==ts,treatment])
    dataset[dataset[,'TS']==ts,treatment] <- (dataset[dataset[,'TS']==ts,treatment] - mode) / std
}
#dataset <- dup_data[!is.na(dup_data[, treatment]),]
#dataset[dataset[,'TS']==3,treatment] <- dataset[dataset[,'TS']==3,treatment] + 100
dataset <- read.csv('simulate_data.csv', stringsAsFactors=F)
matrix_list <- toy_prob(dataset, j_class, 1, treatment)
D <- matrix_list$D
dvec <- matrix_list$dvec
A <- matrix_list$A
b_0 <- matrix_list$b_0
unique_seqs <- matrix_list$unique_seqs
sum_y02 <- matrix_list$sum_y02
sol <- solve.QP(D, dvec, A, b_0)
solution <- sol$solution
sigma_vec <- c(1, 1 / solution[1:4])
mu_vec <- c(0, solution[5:8] / solution[1:4])
# true vec
#true_vec <- c(1/true_sigma[-1], true_mu[-1]/true_sigma[-1], true_theta)
# print obj value
print ('optim vs true on D')
print (0.5 * t(solution) %*% D %*% solution - t(dvec) %*% solution + sum_y02)
print (true_mu_vec)
print (mu_vec)
print (true_sigma_vec)
print (sigma_vec)
#print (0.5 * t(true_vec) %*% D %*% true_vec - t(dvec) %*% true_vec + sum_y02)


# construct norm result table
j_class <- 1:5
theta_vec <- solution[9:length(solution)]
norm.data <- cbind(dataset[,c('TS', 'spot', 'seq')], matrix(rep(NA, 6 * dim(dataset)[1]), ncol=6))
colnames(norm.data) <- c('TS', 'spot', 'seq', 'unique_seq_idx', 'raw', 'theta', 'eps', 'mu', 'sigma')
for (n in 1:dim(dataset)[1]) {
    idx <- return_idx(norm.data[n, 'seq'], unique_seqs)
    norm.data[n, 'raw'] <- dataset[n, treatment]
    norm.data[n, 'unique_seq_idx'] <- idx
    norm.data[n, 'theta'] <- theta_vec[idx]
    norm.data[n, 'sigma'] <- sigma_vec[which(j_class == norm.data[n, 'TS'])]
    norm.data[n, 'mu'] <- mu_vec[which(j_class == norm.data[n, 'TS'])]
    norm.data[n, 'eps'] <- (norm.data[n, 'raw'] - norm.data[n, 'mu']) / norm.data[n, 'sigma'] - norm.data[n, 'theta']
}
write.csv(norm.data, "toy_norm_result_table.csv", row.names=F)
