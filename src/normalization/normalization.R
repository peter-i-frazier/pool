purify_data <- function(data, cutoff) {
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
    new_data <- data[data < (new_M + cutoff * new_sigma)]
    return (new_data)
}

