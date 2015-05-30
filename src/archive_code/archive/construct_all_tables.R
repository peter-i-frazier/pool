rm(list=ls())
source("normalization_util.R")
# quadratic approach
require(quadprog)
root_path = "../.."
org.data <- read.csv(paste(root_path, "/data/all_TS_raw_reading.csv", sep=''), stringsAsFactors=F)
anchor_data <- read.csv(paste(root_path, "/data/anchor_data.csv", sep=''), stringsAsFactors=F)
nolabel_table <- read.csv(paste(root_path, "/data/nolabel_test.csv", sep=''), stringsAsFactors=F)
positive_label_table <- read.csv(paste(root_path, "/data/positive_label.csv", sep=''), stringsAsFactors=F)
not_unlabel_table <- read.csv(paste(root_path, "/data/not_unlabel_table.csv", sep=''), stringsAsFactors=F)
const <- 1.0
delta <- 1e-4
# for sfp
k_class <- c('sfp_1', 'sfp_2')
j_class <- 1:5
group_name <- 'sfp'
dataset <- org.data
matrix_list <- prepare_matrix(dataset, k_class, j_class, group_name, nolabel_table, not_unlabel_table, positive_label_table, const, delta)
D <- matrix_list$D
A <- matrix_list$A
b_0 <- matrix_list$b_0
unique_seqs <- matrix_list$unique_seqs
num_j <- length(j_class)
num_k <- length(k_class)
eps <- 1e-8
sol <- solve.QP(Dmat=D+eps*diag(dim(D)[1]), dvec=rep(0,dim(D)[1]), Amat=A, bvec=b_0)
theta_vec <- sol$solution[(2 * num_j * num_k + 1):length(sol$solution)]
sigma_vec <- 1 / sol$solution[1:(num_j * num_k)]
mu_vec <- sigma_vec * sol$solution[(num_j * num_k + 1):(2 * num_j * num_k)]
result_sfp <- list(unique_seqs=unique_seqs, theta_vec=theta_vec, sigma_vec=sigma_vec, mu_vec=mu_vec, j_class=j_class, k_class=k_class)

# for AcpS
k_class <- c('AcpS_1', 'AcpS_2')
j_class <- c(1, 3, 4, 5)
group_name <- 'AcpS'
dataset <- org.data[org.data[, 'TS'] != 2, ]
matrix_list <- prepare_matrix(dataset, k_class, j_class, group_name, nolabel_table, not_unlabel_table, positive_label_table, const, delta)
D <- matrix_list$D
A <- matrix_list$A
b_0 <- matrix_list$b_0
unique_seqs <- matrix_list$unique_seqs
num_j <- length(j_class)
num_k <- length(k_class)
eps <- 1e-8
sol <- solve.QP(Dmat=D+eps*diag(dim(D)[1]), dvec=rep(0,dim(D)[1]), Amat=A, bvec=b_0)
theta_vec <- sol$solution[(2 * num_j * num_k + 1):length(sol$solution)]
sigma_vec <- 1 / sol$solution[1:(num_j * num_k)]
mu_vec <- sigma_vec * sol$solution[(num_j * num_k + 1):(2 * num_j * num_k)]
result_AcpS <- list(unique_seqs=unique_seqs, theta_vec=theta_vec, sigma_vec=sigma_vec, mu_vec=mu_vec, j_class=j_class, k_class=k_class)

# construct theta_table
theta_table <- data.frame(result_sfp$unique_seqs)
theta_table <- cbind(theta_table, matrix(result_sfp$theta_vec, ncol=2, byrow=T), matrix(rep(NA, length(result_sfp$theta_vec)), ncol=2))
colnames(theta_table) <- c('unique_seqs', 'sfp_1', 'sfp_2', 'AcpS_1', 'AcpS_2')
for (i in 1:length(result_AcpS$unique_seqs)) {
    seq <- result_AcpS$unique_seqs[i]
    i_idx <- return_idx(seq, result_sfp$unique_seqs)
    if (i_idx == 0) {
        print ("error, seq not found in sfp seqs")
    } else {
        theta_table[i_idx, 'AcpS_1'] <- result_AcpS$theta_vec[2 * (i-1) + 1]
        theta_table[i_idx, 'AcpS_2'] <- result_AcpS$theta_vec[2 * (i-1) + 2]
    }
}
write.csv(theta_table, "theta_table.csv", row.names=F)

# construct norm result table
norm.data <- cbind(org.data[,c('TS', 'spot', 'seq')], matrix(rep(NA, 24 * dim(org.data)[1]), ncol=24))
colnames(norm.data) <- c('TS', 'spot', 'seq', 'unique_seq_idx', 'sfp_1_raw', 'sfp_1_theta', 'sfp_1_eps', 'sfp_1_mu', 'sfp_1_sigma', 'sfp_2_raw', 'sfp_2_theta', 'sfp_2_eps', 'sfp_2_mu', 'sfp_2_sigma', 'AcpS_1_raw', 'AcpS_1_theta', 'AcpS_1_eps', 'AcpS_1_mu', 'AcpS_1_sigma', 'AcpS_2_raw', 'AcpS_2_theta', 'AcpS_2_eps', 'AcpS_2_mu', 'AcpS_2_sigma')
for (n in 1:dim(org.data)[1]) {
    norm.data[n, c('sfp_1_theta', 'sfp_2_theta', 'AcpS_1_theta', 'AcpS_2_theta')] <- theta_table[norm.data[n, 'seq'] == theta_table[,'unique_seqs'], c('sfp_1','sfp_2', 'AcpS_1', 'AcpS_2')]
    norm.data[n, c('sfp_1_raw', 'sfp_2_raw', 'AcpS_1_raw', 'AcpS_2_raw')] <- org.data[n, c('sfp_1','sfp_2', 'AcpS_1', 'AcpS_2')]
    norm.data[n, 'unique_seq_idx'] <- return_idx(norm.data[n, 'seq'], result_sfp$unique_seqs)
    norm.data[n, 'sfp_1_mu'] <- result_sfp$mu_vec[(norm.data[n, 'TS']-1) * 2 + 1]
    norm.data[n, 'sfp_1_sigma'] <- result_sfp$sigma_vec[(norm.data[n, 'TS']-1) * 2 + 1]
    norm.data[n, 'sfp_1_eps'] <- (org.data[n, 'sfp_1'] - norm.data[n, 'sfp_1_mu']) / norm.data[n, 'sfp_1_sigma'] - norm.data[n, 'sfp_1_theta']
    norm.data[n, 'sfp_2_mu'] <- result_sfp$mu_vec[(norm.data[n, 'TS']-1) * 2 + 2]
    norm.data[n, 'sfp_2_sigma'] <- result_sfp$sigma_vec[(norm.data[n, 'TS']-1) * 2 + 2]
    norm.data[n, 'sfp_2_eps'] <- (org.data[n, 'sfp_2'] - norm.data[n, 'sfp_2_mu']) / norm.data[n, 'sfp_2_sigma'] - norm.data[n, 'sfp_2_theta']
    if (norm.data[n, 'TS'] == 1) {
        norm.data[n, 'AcpS_1_mu'] <- result_AcpS$mu_vec[(norm.data[n, 'TS']-1) * 2 + 1]
        norm.data[n, 'AcpS_1_sigma'] <- result_AcpS$sigma_vec[(norm.data[n, 'TS']-1) * 2 + 1]
        norm.data[n, 'AcpS_1_eps'] <- (org.data[n, 'AcpS_1'] - norm.data[n, 'AcpS_1_mu']) / norm.data[n, 'AcpS_1_sigma'] - norm.data[n, 'AcpS_1_theta']
        norm.data[n, 'AcpS_2_mu'] <- result_AcpS$mu_vec[(norm.data[n, 'TS']-1) * 2 + 2]
        norm.data[n, 'AcpS_1_sigma'] <- result_AcpS$sigma_vec[(norm.data[n, 'TS']-1) * 2 + 2]
        norm.data[n, 'AcpS_2_eps'] <- (org.data[n, 'AcpS_2'] - norm.data[n, 'AcpS_2_mu']) / norm.data[n, 'AcpS_1_sigma'] - norm.data[n, 'AcpS_2_theta']
    }
    if (norm.data[n, 'TS'] > 2) {
        norm.data[n, 'AcpS_1_mu'] <- result_AcpS$mu_vec[(norm.data[n, 'TS']-2) * 2 + 1]
        norm.data[n, 'AcpS_1_sigma'] <- result_AcpS$sigma_vec[(norm.data[n, 'TS']-2) * 2 + 1]
        norm.data[n, 'AcpS_1_eps'] <- (org.data[n, 'AcpS_1'] - norm.data[n, 'AcpS_1_mu']) / norm.data[n, 'AcpS_1_sigma'] - norm.data[n, 'AcpS_1_theta']
        norm.data[n, 'AcpS_2_mu'] <- result_AcpS$mu_vec[(norm.data[n, 'TS']-2) * 2 + 2]
        norm.data[n, 'AcpS_1_sigma'] <- result_AcpS$sigma_vec[(norm.data[n, 'TS']-2) * 2 + 2]
        norm.data[n, 'AcpS_2_eps'] <- (org.data[n, 'AcpS_2'] - norm.data[n, 'AcpS_2_mu']) / norm.data[n, 'AcpS_1_sigma'] - norm.data[n, 'AcpS_2_theta']
    }
}
write.csv(norm.data, "norm_result_table.csv", row.names=F)
