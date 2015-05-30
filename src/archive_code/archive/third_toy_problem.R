rm(list=ls())
source("normalization_util.R")
source("simulate_data.R")

pt <- function(ts, spot) { print (org.data[org.data[,'TS']==ts & org.data[,'spot']==spot,])}

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
    len_x <- num_j + num_i
    # construct D matrix and dvec
    D_matrix <- matrix(rep(0, len_x * len_x), nrow=len_x, ncol=len_x)
    dvec <- rep(0, len_x)
    sum_y02 <- 0
    for (n in 1:dim(data)[1]) {
        idx_theta <- num_j + y_index[n, 'i']
        D_matrix[idx_theta, idx_theta] <- 1.0 + D_matrix[idx_theta, idx_theta]
        if (y_index[n, 'j'] > 0) {
            idx_alpha <- y_index[n, 'j']
            D_matrix[idx_alpha, idx_alpha] <- data[n, treatment]^2 + D_matrix[idx_alpha, idx_alpha]
            D_matrix[idx_alpha, idx_theta] <- -data[n, treatment] + D_matrix[idx_alpha, idx_theta]
            D_matrix[idx_theta, idx_alpha] <- -data[n, treatment] + D_matrix[idx_theta, idx_alpha]
        } else {
            dvec[idx_theta] <- dvec[idx_theta] + 2 * data[n, treatment]
            sum_y02 <- sum_y02 + data[n, treatment]^2
        }
    }
    tA <- cbind(diag(num_j), matrix(0, nrow=(num_j), ncol=(num_i)))
    b_0 <- rep(1, num_j)
    return (list(D=2*D_matrix, dvec=dvec, A=t(tA), b_0=b_0, sum_y02=sum_y02, unique_seqs=unique_seqs, meq=num_j))
}

# quadratic approach
require(quadprog)
root_path = "../.."
org.data <- read.csv(paste(root_path, "/data/all_TS_raw_reading.csv", sep=''), stringsAsFactors=F)
background.data <- read.csv('background.csv', stringsAsFactors=F)
for (r in 1:nrow(background.data))
    background.data[r, 'val'] <- org.data[org.data[,'TS'] == background.data[r, 'TS'] & 
                                          org.data[,'spot'] == background.data[r, 'spot'], 
                                          background.data[r, 'treatment']]
mean.table <- c()
for (ts in 1:5)
    for (k in c('sfp_1', 'sfp_2', 'AcpS_1', 'AcpS_2'))
        if (!(ts==2 && (k=='AcpS_1' || k=='AcpS_2')))
            mean.table <- rbind(mean.table, data.frame(TS=ts, treatment=k, 
                                mean=mean(background.data[background.data[,'TS']==ts & background.data[,'treatment']==k, 'val']), 
                                stringsAsFactors=F))

treatment <- "sfp_1"
group <- 1
if (group == 1) {
    j_class <- c(1,2,3,5)
} else {
    j_class <- c(1,3,5)
}
dataset <- org.data[!is.na(org.data[, treatment]),]
spread <- 1.25 * (max(dataset[dataset[,'TS']==3,treatment]) - min(dataset[dataset[,'TS']==3,treatment]))
for (ts in 1:5) {
    #std <- sd(dataset[dataset[,'TS']==ts,treatment])
    mu <- min(dataset[dataset[,'TS']==ts,treatment])
    #mu <- mean.table[mean.table[,'TS'] == ts & mean.table[, 'treatment'] == treatment, 'mean']
    #std <- ifelse(ts == 1, 45683 - mu, max(dataset[dataset[,'TS']==ts,treatment]) - mu)
    std <- spread * ifelse(ts == 2, 1/175.718, 1)
    print (std)
    #std <- 45612 - min(dataset[dataset[,'TS']==1,treatment])
    dataset[dataset[,'TS']==ts,treatment] <- (dataset[dataset[,'TS']==ts,treatment] - mu) / std
}
matrix_list <- toy_prob(dataset, j_class, 4, treatment)
D <- matrix_list$D
dvec <- matrix_list$dvec
A <- matrix_list$A
b_0 <- matrix_list$b_0
unique_seqs <- matrix_list$unique_seqs
sum_y02 <- matrix_list$sum_y02
sol <- solve.QP(D, dvec, A, b_0, meq=matrix_list$meq)
#sol <- solve.QP(D, dvec, A, b_0)
solution <- sol$solution

# construct norm result table
if (group == 1) {
    j_class <- 1:5
    theta_vec <- solution[5:length(solution)]
    sigma_vec <- c(1 / solution[1:3], 1, 1 / solution[4])
} else {
    j_class <- c(1,3,4,5)
    theta_vec <- solution[4:length(solution)]
    sigma_vec <- c(1 / solution[1], 1, 1 / solution[2:3])
}
norm.data <- cbind(dataset[,c('TS', 'spot', 'seq')], matrix(NA, nrow(dataset), 5))
colnames(norm.data) <- c('TS', 'spot', 'seq', 'unique_seq_idx', 'raw', 'theta', 'eps', 'sigma')
for (n in 1:nrow(dataset)) {
    idx <- return_idx(norm.data[n, 'seq'], unique_seqs)
    norm.data[n, 'raw'] <- dataset[n, treatment]
    norm.data[n, 'unique_seq_idx'] <- idx
    norm.data[n, 'theta'] <- theta_vec[idx]
    norm.data[n, 'sigma'] <- sigma_vec[which(j_class == norm.data[n, 'TS'])]
    norm.data[n, 'eps'] <- norm.data[n, 'raw'] / norm.data[n, 'sigma'] - norm.data[n, 'theta']
}
write.csv(norm.data, sprintf("%s_norm_result_table.csv", treatment), row.names=F)

dup_data <- c()
for (seq in unique_seqs) {
    sub_data <- norm.data[norm.data[,'seq'] == seq,]
    if (nrow(sub_data) > 1) {
        dup_data <- rbind(dup_data, sub_data)
    }
}
write.csv(dup_data, sprintf("%s_dup_data.csv", treatment), row.names=F)
