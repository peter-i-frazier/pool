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

prepare_nolabel_idx <- function(data, nolabel_table, k_class, unique_seqs) {
    i_idx <- c()
    k_idx <- c()
    for (k in 1:length(k_class)) {
        local_unique_seqs <- c()
        for (n in 1:dim(nolabel_table)[1]) {
            if (nolabel_table[n, 'type'] == k_class[k]) {
                temp <- data[data[,'TS'] == nolabel_table[n, 'TS'],]
                seq <- temp[temp[,'spot'] == paste(nolabel_table[n, 'row'], as.character(nolabel_table[n, 'col']), sep=''), 'seq']
                if (return_idx(seq, local_unique_seqs) == 0) {
                    local_unique_seqs <- c(local_unique_seqs, seq)
                    i_idx <- c(i_idx, return_idx(seq, unique_seqs))
                    k_idx <- c(k_idx, k)
                }
            }
        }
    }
    return (list(i_idx=i_idx, k_idx=k_idx))
}

prepare_positive_label_idx <- function(data, positive_label_table, k_class, unique_seqs) {
    i_idx <- c()
    k_idx <- c()
    for (k in 1:length(k_class)) {
        local_unique_seqs <- c()
        for (n in 1:dim(positive_label_table)[1]) {
            if (positive_label_table[n, 'type'] == k_class[k]) {
                temp <- data[data[,'TS'] == positive_label_table[n, 'TS'],]
                seq <- temp[temp[,'spot'] == paste(positive_label_table[n, 'row'], as.character(positive_label_table[n, 'col']), sep=''), 'seq']
                if (return_idx(seq, local_unique_seqs) == 0) {
                    local_unique_seqs <- c(local_unique_seqs, seq)
                    i_idx <- c(i_idx, return_idx(seq, unique_seqs))
                    k_idx <- c(k_idx, k)
                }
            }
        }
    }
    return (list(i_idx=i_idx, k_idx=k_idx))
}

prepare_not_unlabel_idx <- function(data, not_unlabel_table, unique_seqs, group_name) {
    i_idx <- c()
    local_unique_seqs <- c()
    for (n in 1:dim(not_unlabel_table)[1]) {
        if (not_unlabel_table[n, 'type'] == group_name) {
            temp <- data[data[,'TS'] == not_unlabel_table[n, 'TS'],]
            seq <- temp[temp[,'spot'] == paste(not_unlabel_table[n, 'row'], as.character(not_unlabel_table[n, 'col']), sep=''), 'seq']
            if (return_idx(seq, local_unique_seqs) == 0) {
                local_unique_seqs <- c(local_unique_seqs, seq)
                i_idx <- c(i_idx, return_idx(seq, unique_seqs))
            }
        }
    }
    return (i_idx)
}

prepare_matrix <- function(data, k_class, j_class, group_name, nolabel_table, not_unlabel_table, positive_label_table, const, delta) {
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
    # # sum of theta = constant
    # A_sum_theta <- c(rep(0, 2 * num_j * num_k), rep(1, dim(D_matrix)[1] - 2 * num_j * num_k))
    # theta for positive_label >= const
    positive_label_idx_list <- prepare_positive_label_idx(data, positive_label_table, k_class, unique_seqs)
    positive_label_i_idx <- positive_label_idx_list$i_idx
    positive_label_k_idx <- positive_label_idx_list$k_idx
    A_positive_label <- matrix(rep(0, length(positive_label_i_idx) * len_x), nrow=length(positive_label_i_idx), ncol=len_x)
    for (n in 1:length(positive_label_i_idx)) {
        A_positive_label[n, 2 * num_j * num_k + (positive_label_i_idx[n] - 1) * num_k + positive_label_k_idx[n]] <- 1
    }
    # theta for nolabel <= \delta
    nolabel_idx_list <- prepare_nolabel_idx(data, nolabel_table, k_class, unique_seqs)
    nolabel_i_idx <- nolabel_idx_list$i_idx
    nolabel_k_idx <- nolabel_idx_list$k_idx
    print (length(nolabel_i_idx))
    A_nolabel <- matrix(rep(0, length(nolabel_i_idx) * len_x), nrow=length(nolabel_i_idx), ncol=len_x)
    for (n in 1:length(nolabel_i_idx)) {
        A_nolabel[n, 2 * num_j * num_k + (nolabel_i_idx[n] - 1) * num_k + nolabel_k_idx[n]] <- -1
    }
    # theta_ik1 - theta_ik2 <= \delta for not_unlabel
    not_unlabel_i_idx <- prepare_not_unlabel_idx(data, not_unlabel_table, unique_seqs, group_name)
    A_not_unlabel <- matrix(rep(0, length(not_unlabel_i_idx) * len_x), nrow=length(not_unlabel_i_idx), ncol=len_x)
    for (n in 1:length(not_unlabel_i_idx)) {
        A_not_unlabel[n, 2 * num_j * num_k + (not_unlabel_i_idx[n] - 1) * num_k + 1] <- -1
        A_not_unlabel[n, 2 * num_j * num_k + (not_unlabel_i_idx[n] - 1) * num_k + 2] <- 1
    }
    # meq <- 1 + length(not_unlabel_i_idx)
    # inequality part
    # First 2 * num_j * num_k entries are positive
    A_positive <- matrix(rep(0, 2 * num_j * num_k * len_x), nrow=2 * num_j * num_k, ncol=len_x)
    for (i in 1:(2 * num_j * num_k)) {
        A_positive[i,i] <- 1
    }
    # theta_ik1 - theta_ik2 >= 0
    A_theta <- matrix(rep(0, length(unique_seqs) * len_x), nrow=length(unique_seqs), ncol=len_x)
    for (i in 1:length(unique_seqs)) {
        A_theta[i, 2 * num_j * num_k + (i - 1) * num_k + 1] <- 1
        A_theta[i * num_j * num_k + (i - 1) * num_k + 2] <- -1
    }
    A_matrix <- rbind(A_positive_label, A_nolabel, A_not_unlabel, A_positive, A_theta)
    # b_0
    #b_0 <- c(num_i , rep(0, dim(A_matrix)[1]-1))
    b_0 <- c(rep(const, dim(A_positive_label)[1]), rep(-delta, dim(A_nolabel)[1]), rep(-delta, dim(A_not_unlabel)[1]), rep(0, dim(A_positive)[1]), rep(0, dim(A_theta)[1]))
    # # construct neq matrix: neq_matrix %*% x <= 0
    # neq_matrix <- matrix(rep(0, (2*num_j*num_k+length(unique_seqs))* len_x), nrow=2*num_j*num_k+length(unique_seqs), ncol=len_x)
    # for (i in 1:(2*num_j*num_k)) {
    #     neq_matrix[i,i] <- -1
    # }
    # for (i in 1:length(unique_seqs)) {
    #     neq_matrix[2*num_j*num_k+i, 2 * num_j * num_k + (i - 1) * num_k + 1] <- -1
    #     neq_matrix[2*num_j*num_k+i, 2 * num_j * num_k + (i - 1) * num_k + 2] <- 1
    # }
    # # construct eq matrix: eq_matrix %*% x = 0
    # eq_matrix <- A_unlabel
    return (list(D=D_matrix, A=t(A_matrix), b_0=b_0, unique_seqs=unique_seqs))
}


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
