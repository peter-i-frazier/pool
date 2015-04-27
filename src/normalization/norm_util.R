require(quadprog)

pt <- function(ts, spot) { print (org.data[org.data[,'TS']==ts & org.data[,'spot']==spot,])}

return.idx <- function(seq, seq.vec) {
    if (length(seq.vec) > 0)
        for (i in 1:length(seq.vec))
            if (seq == seq.vec[i])
                return (i)
    return (0)
}

prepare.matrix <- function(data, const, treatment, j.class) {
# Construct matrix D and A for quadratic program with linear constraint, which has form
# min -d^T b + 1/2 x^T D x, s.t A^T x >= b.0. Note we substitute constraint \sum \theta_i = const
# into objective by \theta_N = const - \sum_{i=1:N-1} \theta_i
# Input: 
#       data: data for duplicated peptides
#       const:\sum \theta_i = const
#       treatment: choose from "sfp_1", "sfp_2", "AcpS_1" or "AcpS_2"
# Output:
#       D, d, A, b.0, dup_unique_seqs
    ij_idx <- matrix(NA, nrow=nrow(data), ncol=2)
    colnames(ij_idx) <- c('i', 'j')
    unique_seqs <- c()
    for (n in 1:nrow(data)) {
        if (return.idx(data[n, 'seq'], unique_seqs) == 0)
            unique_seqs <- c(unique_seqs, data[n, 'seq'])
        ij_idx[n, 'i'] <- return.idx(data[n, 'seq'], unique_seqs)
        ij_idx[n, 'j'] <- data[n, 'TS']
    }
    N <- length(unique_seqs)
    len_x <- length(j.class) + N - 1
    # construct D matrix and dvec
    D_matrix <- matrix(0, nrow=len_x, ncol=len_x)
    dvec <- rep(0, len_x)
    for (n in 1:nrow(data)) {
        beta_idx <- which(j.class == ij_idx[n, 'j'])
        theta_idx <- ij_idx[n, 'i'] + length(j.class)
        y <- data[n, treatment]
        D_matrix[beta_idx, beta_idx] <- D_matrix[beta_idx, beta_idx] + y^2
        if (ij_idx[n, 'i'] != N) {
            D_matrix[beta_idx, theta_idx] <- D_matrix[beta_idx, theta_idx] - y
            D_matrix[theta_idx, beta_idx] <- D_matrix[theta_idx, beta_idx] - y
            D_matrix[theta_idx, theta_idx] <- D_matrix[theta_idx, theta_idx] + 1
        } else {
            for (i in 1:(N-1)) {
                for (j in 1:(N-1)) {
                    D_matrix[length(j.class) + i, length(j.class) + j] <- D_matrix[length(j.class) + i, length(j.class) + j] + 0.5
                    D_matrix[length(j.class) + j, length(j.class) + i] <- D_matrix[length(j.class) + j, length(j.class) + i] + 0.5
                }
                dvec[length(j.class) + i] <- dvec[length(j.class) + i] + 2 * const
                D_matrix[beta_idx, length(j.class) + i] <- D_matrix[beta_idx, length(j.class) + i] + y
                D_matrix[length(j.class) + i, beta_idx] <- D_matrix[length(j.class) + i, beta_idx] + y
            }
            dvec[beta_idx] <- dvec[beta_idx] + 2 * const * y
        }
    }
    tA <- cbind(diag(length(j.class)), matrix(0, nrow=length(j.class), ncol=(N-1)))
    b.0 <- rep(0, length(j.class))
    return (list(D=2*D_matrix, d=dvec, A=t(tA), b.0=b.0, dup_unique_seqs=unique_seqs))
}

normalize <- function(org.data, dup.data, treatment, result.path=NA, const=NA) {
    if (is.na(const))
        const <- length(unique(dup.data[, 'seq'])) - 1
    matrix.list <- prepare.matrix(dup.data, const, treatment, j.class)
    D <- matrix.list$D
    d <- matrix.list$d
    A <- matrix.list$A
    b.0 <- matrix.list$b.0
    dup.unique.seqs <- matrix.list$dup_unique_seqs
    sol <- solve.QP(D, d, A, b.0) 
    solution <- sol$solution
# end solve QP
    if (!is.na(result.path)) {
# construct norm table
        N <- length(dup.unique.seqs)
        norm.data <- cbind(org.data[,c('TS', 'spot', 'seq')], matrix(NA, nrow(org.data), 5))
        colnames(norm.data) <- c('TS', 'spot', 'seq', 'unique_seq_idx', 'raw', 'theta', 'eps', 'sigma')
        for (n in 1:nrow(org.data)) {
            norm.data[n, 'unique_seq_idx'] <- return.idx(norm.data[n, 'seq'], unique.seqs)
            norm.data[n, 'raw'] <- org.data[n, treatment]
            norm.data[n, 'sigma'] <- 1 / solution[which(j.class == norm.data[n, 'TS'])]
            dup.idx <- return.idx(norm.data[n, 'seq'], dup.unique.seqs)
            if (dup.idx == 0) {
                norm.data[n, 'theta'] <- norm.data[n, 'raw'] / norm.data[n, 'sigma']
            } else if (dup.idx == N) {
                norm.data[n, 'theta'] <- const - sum(solution[(length(j.class)+1):length(solution)])
            } else {
                norm.data[n, 'theta'] <- solution[length(j.class) + dup.idx]
            }
            norm.data[n, 'eps'] <- norm.data[n, 'raw'] / norm.data[n, 'sigma'] - norm.data[n, 'theta']
        }
        multiplier <- 1 / max(norm.data[, 'theta'])
        norm.data[, c('theta', 'eps')] <- norm.data[, c('theta', 'eps')] * multiplier
        norm.data[, 'sigma'] <- norm.data[, 'sigma'] / multiplier
        write.csv(norm.data, sprintf("%s/%s_norm_table.csv", result.path, treatment), row.names=F)
# construct dup table
        dup.data <- c()
        for (seq in dup.unique.seqs) {
            dup.data <- rbind(dup.data, norm.data[norm.data[, 'seq'] == seq,])
        }
        write.csv(dup.data, sprintf("%s/%s_dup_table.csv", result.path, treatment), row.names=F)
    }
    return (solution)
}
