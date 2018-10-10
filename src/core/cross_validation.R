
leave_n_out_cv <- function(dataset_X, dataset_Y, alpha.1, alpha.0, p1, mc_itr, num_out) {
  # N-fold cross validation
  #
  # Args:
  #   dataset_X: matrix, each row is feature vector for one peptide
  #   dataset_Y: vector, 0 or 1
  #   alpha.1: list of vectors, prior
  #   alpha.0: list of vectors, prior
  #   p1: double, P(Y=1)
  #   mc_itr: int, number of MC iterations for estimating probability
  #   num_out: int, number to leave out in cross validation
  #
  # Returns a list:
  #   prob:
  #   Y:
    num_data <- dim(dataset_X)[1]
    if (num_data <= 1 || num_out > num_data)
      stop("ERROR: number of data <= 1 or num_out > num_data")
    num_groups <- floor(num_data / num_out)
    prob <- c()
    for (n in 1:(num_groups + 1)) {
      if (n <= num_groups) {
        train_X <- dataset_X[-(((n-1) * num_out + 1):(n * num_out)),]
        train_Y <- dataset_Y[-(((n-1) * num_out + 1):(n * num_out))]
        test_X <- dataset_X[(((n-1) * num_out + 1):(n * num_out)),]
      } else if (num_groups * num_out != num_data) {
        train_X <- dataset_X[-(((n-1) * num_out + 1):num_data),]
        train_Y <- dataset_Y[-(((n - 1) * num_out + 1):num_data)]
        test_X <- dataset_X[(((n - 1) * num_out + 1):num_data),]
      } else {
        break
      }
      if (is.null(nrow(test_X))) {
        test_X <- matrix(test_X, nrow=1, ncol=length(test_X))
      }
      prob <- c(prob, Predict(train_X, train_Y, test_X, alpha.1, alpha.0, p1, mc_itr))
    }
    return (list(prob=prob, Y=dataset_Y))
}

GetRocXy <- function(prob, Y) {
  # Prepare x,y values for ROC plot
  # Returns data.frame
  N <- length(prob)
  FPR <- rep(-1, N)
  TPR <- rep(-1, N)
  thresholds <- sort(prob)
  for( i in 1:N) {
    threshold <- thresholds[i]
    label <- rep(0, N)
    for( j in 1:N){
      if(prob[j] >= threshold) {
        label[j] <- 1 }
    }
    FPR[N - i + 1] <- sum((Y==0)&(label==1))/sum(Y==0)
    TPR[N - i + 1] <- sum((Y==1)&(label==1))/sum(Y==1)
  }
  return (data.frame(x = FPR, y = TPR))
}

GetAuc<- function(x, y) {
# Calculate area under curve
  U <- 0
  L <- 0
  for(i in 1:(length(y) - 1)) {
    U <- U + y[i+1] * (x[i+1] - x[i])
    L <- L + y[i] * (x[i+1] - x[i])
  }
  return ((U + L) / 2)
}

