# Normalize the raw data
# Inputs: raw data (vector)
# Outputs: normalized data (vector)
normalize_data <- function(data, cutoff) {
    M <- median(data)
    L <- quantile(data, 0.25)
    z <- qnorm(0.25)
    sigma <- (M - L) / (-z)
    processed_data <- c()
    for (i in 1:length(data)) {
        if (data[i] < (M + cutoff*sigma)) {
            processed_data <- c(processed_data, data[i])
        }
    }
    new_M <- median(processed_data)
    new_L <- quantile(processed_data, 0.25)
    new_sigma <- (new_M - new_L) / (-z)
    ndata <- (data - new_M) / new_sigma
    print ('new M')
    print (new_M)
    print ('new sigma')
    print (new_sigma)
    return (ndata)
}

normalize_data_manual <- function(data, M) {
    L <- quantile(data, 0.05)
    z <- qnorm(0.05)
    sigma <- (M - L) / (-z)
    ndata <- (data - M) / sigma
    print ('M')
    print (M)
    print ('sigma')
    print (sigma)
    return (ndata)
}

# Create indicator vectors for peptides
# Inputs: 
#       normalized data: N * 4 matrix (N is no. of peptides, 4 columns are: 
#             1. Intensity after sfp treatment; 
#             2. Intensity after sfp_PfAcpH;
#             3. Intensity after AcpS; 
#             4. Intensity after AcpS_PfAcpH
#       threshold: vector with length 4; each element is threshold for corresponding column of data
# Output: 
#       N * 4 matrix (each row is indicator vector for nth peptide)
peptide_indicator <- function(normal_data, threshold) {
    to_return <- c()
    for (i in 1:4) {
        to_return <- cbind(to_return, (normal_data[,i] > threshold[i]))
    }
    to_return <- 1 * to_return
    return (to_return)
}

# Put peptides in categories, specified in README
# Input:
#       indicator: N * 4 matrix, which is output from function peptide_indicator(...)
# Output:
#       vector with length N, tell which category is nth peptide in (1-7), none of the category if -1
categorize_peptide <- function(indicator) {
    N <- dim(indicator)[1]
    table <- matrix(c(0,1,1,1,1,0,0,0,0,1,0,1,0,0,0,1,1,0,0,1,1,0,0,1,0,0,0,1), nrow = 7)
    to_return <- rep(-1, N)
    for (n in 1:N) {
        for (i in 1:7) {
            if (prod(indicator[n,] == table[i,]) == 1) {
                to_return[n] = i
                break
            }
        }
    }
    return (to_return)
}

# convert code_idx -> spot_idx  &  spot_idx -> code_idx
# spot_idx is the position of a peptide in the spot array, for example, "B10"
# code_idx is index of peptide we used in code, we treat spot_idx as a row-major matrix, so "B10"= (2-1)*30 + 10 = 40
# Input:
#       vector, and each element stores an index
# Output:
#       vector, each element stores corresponding converted index
convert_to_spot_idx <- function(index) {
    row_letter <- c('A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S')
    spot_idx <- c()
    for (i in 1:length(index)) {
        if (index[i] %% 30 == 0) {
            spot_idx <- c(spot_idx, paste(row_letter[index[i]%/%30], 30, sep=''))
        }
        else {
            spot_idx <- c(spot_idx, paste(row_letter[index[i]%/%30 + 1], index[i]%%30, sep=''))
        }
    }
    return (spot_idx)
}

convert_to_code_idx <- function(index) {
    row_letter <- c('A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S')
    N <- length(index)
    code_idx <- rep(-1,N)
    for (i in 1:N) {
        split_str <- unlist(strsplit(index[i],''))
        a <- which(row_letter == split_str[1])
        b <- as.numeric(paste(split_str[-1], collapse=''))
        code_idx[i] <- (a-1)*30 + b
    }
    return (code_idx)
}



