root_path = "../.."
org.data <- read.csv(paste(root_path, "/data/all_TS_raw_reading.csv", sep=''), stringsAsFactors=F)
anchor_data <- read.csv(paste(root_path, "/data/anchor_data.csv", sep=''), stringsAsFactors=F)
unlabel_table <- read.csv(paste(root_path, "/data/unlabel_table.csv", sep=''), stringsAsFactors=F)

return_idx <- function(seq, seq_list) {
    idx <- 0
    if (length(seq_list) > 0) {
        for (i in 1:length(seq_list)) {
            if (seq == seq_list[i]) {
                idx <- i
                break
            }
        }
    }
    return (idx)
}

prepare_anchor_data <- function(data) {
    anchor_data <- c()
    for (i in 1:dim(data)[1]) {
        for (j in 1:dim(data)[1]) {
            if (data[i, 'seq'] == data[j, 'seq'] && i != j) {
                anchor_data <- rbind(anchor_data, data[i,])
                break
            }
        }
    }
    return (anchor_data)
}

prepare_unlabel_list <- function(data, unlabel_table, k_class) {
    to_return <- list()
    for (k in 1:length(k_class)) {
        unique_seqs <- c()
        for (n in 1:dim(unlabel_table)[1]) {
            if (unlabel_table[n, k_class[k]] != -1) {
                TS_spot <- strsplit(unlabel_table[n, k_class[k]], split='_')[[1]]
                temp <- data[data[,'TS']==as.numeric(TS_spot[1]),]
                seq <- temp[temp[,'spot']==TS_spot[2], 'seq']
                if (return_idx(seq, unique_seqs) == 0) {
                    unique_seqs <- c(unique_seqs, seq)
                }
            }
        }
        to_return[[k_class[k]]] <- unique_seqs
    }
    return (to_return)
}

prepare_matrix <- function(data, unlabel_list, k_class, j_class) {
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
        y_index[n, 'j'] <- data[n, 'TS']
    }
    num_i <- length(unique_seqs)
    num_j <- length(j_class)
    num_k <- length(k_class)
    len_x <- num_j * num_k + num_j * num_k + num_i * num_k
    # construct D matrix
    D_matrix <- matrix(rep(0, len_x * len_x), nrow=len_x, ncol=len_x)
    for (n in 1:dim(data)[1]) {
        for (k in 1:num_k) {
            idx_alpha <- (y_index[n, 'j'] - 1) * num_k + k
            idx_beta <- num_j * num_k + (y_index[n, 'j'] - 1) * num_k + k
            idx_theta <- 2 * num_j * num_k + (y_index[n, 'i'] - 1) * num_k + k
            D_matrix[idx_alpha, idx_alpha] <- data[n, k_class[k]]^2 + D_matrix[idx_alpha, idx_alpha]
            D_matrix[idx_beta, idx_beta] <- 1.0 + D_matrix[idx_beta, idx_beta]
            D_matrix[idx_theta, idx_theta] <- 1.0 + D_matrix[idx_theta, idx_theta]
            D_matrix[idx_alpha, idx_beta] <- -data[n, k_class[k]] + D_matrix[idx_alpha, idx_beta]
            D_matrix[idx_beta, idx_alpha] <- -data[n, k_class[k]] + D_matrix[idx_beta, idx_alpha]
            D_matrix[idx_alpha, idx_theta] <- -data[n, k_class[k]] + D_matrix[idx_alpha, idx_theta]
            D_matrix[idx_theta, idx_alpha] <- -data[n, k_class[k]] + D_matrix[idx_theta, idx_alpha]
            D_matrix[idx_beta, idx_theta] <- 1.0 + D_matrix[idx_beta, idx_theta]
            D_matrix[idx_theta, idx_beta] <- 1.0 + D_matrix[idx_theta, idx_beta]
        }
    }
    # construct A matrix: A_matrix %*% x >= b_0, where first meq are equality constraints
    # equality part
    A_sum_theta <- c(rep(0, 2 * num_j * num_k), rep(1, dim(D)[1] - 2 * num_j * num_k))
    A_unlabel <- c()
    for (k in 1:num_k) {
        unlabel_seqs <- unlabel_list[[k_class[k]]]
        eq_part <- matrix(rep(0, length(unlabel_seqs) * len_x), nrow=length(unlabel_seqs), ncol=len_x)
        for (n in 1:length(unlabel_seqs)) {
            i <- return_idx(unlabel_seqs[n], unique_seqs)
            eq_part[n, 2 * num_j * num_k + (i - 1) * num_k + k] <- 1
        }
        A_unlabel <- rbind(A_unlabel, eq_part)
    }
    meq <- 1 + dim(A_unlabel)[1]
    # inequality part
    A_positive <- matrix(rep(0, 2 * num_j * num_k * len_x), nrow=2 * num_j * num_k, ncol=len_x)
    for (i in 1:(2 * num_j * num_k)) {
        A_positive[i,i] <- 1
    }
    A_theta <- matrix(rep(0, length(unique_seqs) * len_x), nrow=length(unique_seqs), ncol=len_x)
    for (i in 1:length(unique_seqs)) {
        A_theta[i, 2 * num_j * num_k + (i - 1) * num_k + 1] <- 1
        A_theta[i * num_j * num_k + (i - 1) * num_k + 2] <- -1
    }
    A_matrix <- rbind(A_sum_theta, A_unlabel, A_positive, A_theta)
    # b_0
    b_0 <- c(1000, rep(0, dim(A_matrix)[1]-1))
    # construct neq matrix: neq_matrix %*% x <= 0
    neq_matrix <- matrix(rep(0, (2*num_j*num_k+length(unique_seqs))* len_x), nrow=2*num_j*num_k+length(unique_seqs), ncol=len_x)
    for (i in 1:(2*num_j*num_k)) {
        neq_matrix[i,i] <- -1
    }
    for (i in 1:length(unique_seqs)) {
        neq_matrix[2*num_j*num_k+i, 2 * num_j * num_k + (i - 1) * num_k + 1] <- -1
        neq_matrix[2*num_j*num_k+i, 2 * num_j * num_k + (i - 1) * num_k + 2] <- 1
    }
    # construct eq matrix: eq_matrix %*% x = 0
    eq_matrix <- A_unlabel
    return (list(D=D_matrix, A=t(A_matrix), b_0=b_0, neq=neq_matrix, eq=eq_matrix, g=rbind(neq_matrix, eq_matrix, -eq_matrix)))
}


k_class <- c('sfp_step_1', 'sfp_step_2')
j_class <- 1:5
#k_class <- c('AcpS_step_1', 'AcpS_step_2')
#j_class <- c(1, 3, 4, 5)

# # BFGS approach
# require(nloptr)
# unlabel_list <- prepare_unlabel_list(anchor_data, unlabel_table, k_class)
# ll <- prepare_matrix(anchor_data, unlabel_list, k_class, j_class)
# D <- ll$D
# neq_mx <- ll$neq
# eq_mx <- ll$eq
# g_mx <- ll$g
# num_j <- length(j_class)
# num_k <- length(k_class)
# x_0 <- c(rep(30, num_j*num_k), rep(100, num_j*num_k),rep(0, dim(D)[1]-2*num_j*num_k))
# obj_f <- function(x) {
#     xx <- x
#     xx[1:(num_j*num_k)] <- 1 / x[1:(num_j*num_k)]
#     xx[(num_j*num_k+1):(2*num_j*num_k)] = x[(num_j*num_k+1):(2*num_j*num_k)] / x[1:(num_j*num_k)]
#     return (t(xx) %*% D %*% xx)
# }
# print ("start nloptr")
# local_opts = list("algorithm"="NLOPT_LD_LBFGS", "xtol_rel"=1.0e-9)
# #result <- nloptr(x_0, eval_f= (function (x) t(x) %*% D %*% x), eval_grad_f= (function (x) 2 * D %*% x), eval_g_ineq= (function (x) neq_mx %*% x), eval_jac_g_ineq= (function (x) neq_mx), eval_g_eq= (function (x) eq_mx %*% x), eval_jac_g_eq= (function (x) eq_mx), opts=list("algorithm"="NLOPT_LD_AUGLAG", "xtol_rel"=1.0-9, "local_opts"=local_opts))
# result <- nloptr(x_0, eval_f= obj_f, eval_g_ineq= (function (x) g_mx %*% x), eval_jac_g_ineq= (function (x) g_mx), opts=list("algorithm"="NLOPT_LN_COBYLA", "xtol_rel"=1.0-9, "maxeval"=100000))

calculate_norm <- function(data, mu, sigma) {
    norm <- data[,2:6]
    for (i in 1:dim(data)[1]) {
        for (k in 1:length(k_class)) {
            norm[i, k_class[k]] <- (data[i, k_class[k]] - mu[(data[i,'TS']-1) * length(k_class) + k]) / sigma[(data[i,'TS']-1) * length(k_class) + k]
        }
    }
    return (norm)
}

check_peptide <- function(data, seq) {
    result <- c()
    for (i in 1:dim(data)[1]) {
        if (data[i,'seq']==seq) {
            result <- rbind(result, data[i,])
        }
    }
    return (result)
}

sigma <- result$solution[1:10]
mu <- result$solution[11:20]
norm_data <- calculate_norm(anchor_data, mu, sigma)

# quadratic approach
require(quadprog)
unlabel_list <- prepare_unlabel_list(org.data, unlabel_table, k_class)
ll <- prepare_matrix(org.data, unlabel_list, k_class, j_class)
D <- ll$D
A <- ll$A
b_0 <- ll$b_0
num_j <- length(j_class)
num_k <- length(k_class)
eps <- 1e-8
result <- solve.QP(Dmat=D+eps*diag(dim(D)[1]), dvec=rep(0,dim(D)[1]), Amat=A, bvec=b_0)

