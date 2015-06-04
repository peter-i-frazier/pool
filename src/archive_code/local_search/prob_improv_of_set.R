# Compute probability improvement of a set
root_path <- ""
source(paste(root_path,"NB_utility.R",sep=''))

ComputeProbImprovOfSet <- function(S,  # set of recommendations
                                   trainX, trainY, classlist,
                                   S.Pos, maxL, maxR,
                                   gamma.0, gamma.1, prior, iter) {
  print("Computing PI(S)")
  sum <- 0  # cumulates probabilities to get prob improv
  grow.product <- 1
  S.len <- dim(S)[1]
  X <- trainX
  Y <- trainY

  for (i in 1:S.len) {
    test.x <- S[i, ]
    # ptm.0 <- proc.time()

    prob.positive <- Naive_Bayes(X, Y, test.x, classlist, S.Pos,
                                 maxL, maxR, gamma.0, gamma.1, prior, iter)

    sum <- sum + grow.product * prob.positive
    X <- rbind(X, test.x)
    Y <- c(Y, 0)
    grow.product <- grow.product * (1 - prob.positive)

    # print(sprintf("Updating values costs %f", proc.time() - ptm.1))
  }
  return(sum)
}
