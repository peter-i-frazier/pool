rm(list=ls())
source("normalization_util.R")
# quadratic approach
require(quadprog)
root_path = "../.."
org.data <- read.csv(paste(root_path, "/data/all_TS_raw_reading.csv", sep=''), stringsAsFactors=F)
anchor_data <- read.csv(paste(root_path, "/data/anchor_data.csv", sep=''), stringsAsFactors=F)
nolabel_table <- read.csv(paste(root_path, "/data/nolabel_table.csv", sep=''), stringsAsFactors=F)
not_unlabel_table <- read.csv(paste(root_path, "/data/not_unlabel_table.csv", sep=''), stringsAsFactors=F)
const <- 1000
# for sfp
k_class <- c('sfp_1', 'sfp_2')
j_class <- 1:5
group_name <- 'sfp'
dataset <- org.data
matrix_list <- prepare_matrix(dataset, k_class, j_class, group_name, nolabel_table, not_unlabel_table, const)
D <- matrix_list$D
A <- matrix_list$A
b_0 <- matrix_list$b_0
meq <- matrix_list$meq
unique_seqs <- matrix_list$unique_seqs
num_j <- length(j_class)
num_k <- length(k_class)
eps <- 1e-8
sol <- solve.QP(Dmat=D+eps*diag(dim(D)[1]), dvec=rep(0,dim(D)[1]), Amat=A, bvec=b_0, meq=meq)
theta_vec <- sol$solution[(2 * num_j * num_k + 1):length(sol$solution)]
result_sfp <- list(unique_seqs=unique_seqs, theta_vec=theta_vec, j_class=j_class, k_class=k_class)

# for AcpS
k_class <- c('AcpS_1', 'AcpS_2')
j_class <- c(1, 3, 4, 5)
group_name <- 'AcpS'
dataset <- org.data[org.data[, 'TS'] != 2, ]
matrix_list <- prepare_matrix(dataset, k_class, j_class, group_name, nolabel_table, not_unlabel_table, const)
D <- matrix_list$D
A <- matrix_list$A
b_0 <- matrix_list$b_0
meq <- matrix_list$meq
unique_seqs <- matrix_list$unique_seqs
num_j <- length(j_class)
num_k <- length(k_class)
eps <- 1e-8
sol <- solve.QP(Dmat=D+eps*diag(dim(D)[1]), dvec=rep(0,dim(D)[1]), Amat=A, bvec=b_0, meq=meq)
theta_vec <- sol$solution[(2 * num_j * num_k + 1):length(sol$solution)]
result_AcpS <- list(unique_seqs=unique_seqs, theta_vec=theta_vec, j_class=j_class, k_class=k_class)

# construct theta_table
table <- cbind(result_sfp$unique_seqs, matrix(result_sfp$theta_vec, ncol=2, byrow=T), matrix(rep(NA, length(result_sfp$theta_vec)), ncol=2))
colnames(table) <- c('unique_seqs', 'sfp_1', 'sfp_2', 'AcpS_1', 'AcpS_2')
for (i in 1:length(result_AcpS$unique_seqs)) {
    seq <- result_AcpS$unique_seqs[i]
    i_idx <- return_idx(seq, result_sfp$unique_seqs)
    if (i_idx == 0) {
        print ("error, seq not found in sfp seqs")
    } else {
        table[i_idx, 'AcpS_1'] <- result_AcpS$theta_vec[2 * (i-1) + 1]
        table[i_idx, 'AcpS_2'] <- result_AcpS$theta_vec[2 * (i-1) + 2]
    }
}
write.csv(table, "theta_table.csv", row.names=F)
