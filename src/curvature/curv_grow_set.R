# Author: Tom Fei
# Created: 06.07.2015
# Computes curvature for a growing recommendation set

rm(list = ls())
set.seed(1)

root <- "~/repos/peptide-catalysis/src/"
source(paste0(root, "core/dependency.R"))
ResolveDependency(root)

dataset <- GetDataReducedAA("SELECT * FROM binary_labeling_activity")
X <- dataset$feature
Y <- dataset$data[, 'type1']
train.X <- X[Y != -1, ]
train.Y <- Y[Y != -1]

# ====================================================================
# Generate recommendation set
print("Begin generating recommendation set")
gamma.1 <- gamma_1_type1
gamma.0 <- gamma_0_type1
prior.prob <- prior_type1
prior.alpha.1 <- SetPriorReducedAA(gamma.1, NUM_CLASS)
prior.alpha.0 <- SetPriorReducedAA(gamma.0, NUM_CLASS)
NUM.RECOM <- 71
NB.ITER <- 10000
S <- GenRecomSetOld(train.X, train.Y, prior.alpha.1, prior.alpha.0,
                    prior.prob, NUM.RECOM, NB.ITER, minL, maxL, minR, maxR)
print("Done")

# ====================================================================
# Compute curvature
print("Begin computing curvature")
grow.product <- 1
curv.vec <- rep(-1, nrow(S))
X <- train.X
Y <- train.Y
num.neg.val <- 0
num.curv.computed <- 0
while (num.curv.computed < nrow(S)) {
# for (i in 1:nrow(S)) {
  i <- num.curv.computed + 1
  test.x <- S[i, ]
  # ptm.0 <- proc.time()

  alphas <- BayesianNaiveBayes(X, Y, prior.alpha.1, prior.alpha.0, 
                               prior.prob)
  thetas.1 <- SampleThetas(alphas$post.alpha.1, NB.ITER)
  thetas.0 <- SampleThetas(alphas$post.alpha.0, NB.ITER)

  # Note that test.x is a vector. Need to convert it to
  # "row" matrix to pass into CalculateProb.
  prob.list <- CalculateProb(t(test.x), thetas.1, thetas.0, prior.prob)
  prob.positive <- prob.list$mean
  numerator <- grow.product * prob.positive
  denominator <- ComputeProbImprovOfSet(t(test.x), train.X, train.Y,
                                        prior.alpha.1, prior.alpha.0,
                                        prior.prob, NB.ITER)
  curv.vec[i] <- 1 - numerator / denominator
  if (i == 1) 
    curv.vec[i] <- 0
  if (curv.vec[i] < 0) {
    print("Found a negative value!")
    num.neg.val <- num.neg.val + 1
    next
  }
  num.curv.computed <- num.curv.computed + 1
  # Update variables
  X <- rbind(X, test.x)
  Y <- c(Y, 0)
  grow.product <- grow.product * (1 - prob.positive)

  print(sprintf("Curvature for iter #%d is %.3f", i, curv.vec[i]))

}
