# Author: Tom Fei
# Created: 06.08.2015
# Computes curvature for a fixed recommendation set

# ====================================================================
# Initialization
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
class.vec <- QuerySql("sfp_AcpS", "SELECT * FROM reduced_AA")
class.vec <- unique(class.vec[, 'class'])

# ====================================================================
# Generate recommendation set
print("Begin generating recommendation set")
gamma.1 <- gamma_1_type1
gamma.0 <- gamma_0_type1
prior.prob <- prior_type1
prior.alpha.1 <- SetPriorReducedAA(gamma.1, NUM_CLASS)
prior.alpha.0 <- SetPriorReducedAA(gamma.0, NUM_CLASS)
NUM.RECOM <- 71
NB.ITER <- 1000
S <- GenRecomSetOld(train.X, train.Y, prior.alpha.1, prior.alpha.0,
                    prior.prob, NUM.RECOM, NB.ITER, minL, maxL, minR, maxR)
print("Done")

# ====================================================================
# Compute grow product for S
# This is a preliminary step for computing the curvature
print("Begin computing grow product")
grow.product <- 1
X <- train.X
Y <- train.Y

for (i in 1:nrow(S)) {
  test.x <- S[i, ]
  alphas <- BayesianNaiveBayes(X, Y, prior.alpha.1, prior.alpha.0, 
                               prior.prob)
  thetas.1 <- SampleThetas(alphas$post.alpha.1, NB.ITER)
  thetas.0 <- SampleThetas(alphas$post.alpha.0, NB.ITER)

  # Note that test.x is a vector. Need to convert it to
  # "row" matrix to pass into CalculateProb.
  prob.list <- CalculateProb(t(test.x), thetas.1, thetas.0, prior.prob)
  prob.positive <- prob.list$mean
  grow.product <- grow.product * (1 - prob.positive)
  X <- rbind(X, test.x)
  Y <- c(Y, 0)
}
print(sprintf("Grow product for S is %.3f", grow.product))
 
# ====================================================================
# Compute curvature
# Note that now X == rbind(train.X, S)
# and Y == c(train.Y, rep(0, nrow(S)))
# From now on alphas and thetas are fixed
alphas <- BayesianNaiveBayes(X, Y, prior.alpha.1, prior.alpha.0, 
                             prior.prob)
thetas.1 <- SampleThetas(alphas$post.alpha.1, NB.ITER)
thetas.0 <- SampleThetas(alphas$post.alpha.0, NB.ITER)
CV.ITER <- 100
curv.vec <- rep(-1, CV.ITER)
num.neg.val <- 0
num.curv.computed <- 0

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Debug variables
n.vec <- curv.vec
d.vec <- curv.vec
p.vec <- curv.vec
D <- c()
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
print(sprintf("CV.ITER = %d", CV.ITER))

while (num.curv.computed < CV.ITER) {
# for (i in 1:CV.ITER) {
  test.x <- GenOneRandomPeptide(minL, maxL, minR, maxR,
                                MAXL, MAXR, class.vec)

  D <- rbind(D, test.x)
  ptm.0 <- proc.time()

  prob.list <- CalculateProb(t(test.x), thetas.1, thetas.0, prior.prob)
  prob.positive <- prob.list$mean
  numerator <- grow.product * prob.positive
  denominator <- ComputeProbImprovOfSet(t(test.x), train.X, train.Y,
                                        prior.alpha.1, prior.alpha.0,
                                        prior.prob, NB.ITER)
  if (1 - numerator / denominator < 0) {
    print("Found a negative value!")
    num.neg.val <- num.neg.val + 1
    next
  }
  num.curv.computed <- num.curv.computed + 1
  i <- num.curv.computed
  curv.vec[i] <- 1 - numerator / denominator
  p.vec[i] <- prob.positive
  n.vec[i] <- numerator
  d.vec[i] <- denominator

  print(sprintf("Computing curvature costs %1.2f", (proc.time() - ptm.0)[1]))
  print(sprintf("Curvature for iter #%d is %.3f", i, curv.vec[i]))
  print("==============")
}
