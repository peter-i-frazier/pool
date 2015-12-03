# Author: Tom Fei
# Created: 06.06.2015
# Local search algorithm, version 2

rm(list = ls())
set.seed(1)

root <- "~/repos/peptide-catalysis/src/"
source(paste0(root, "core/dependency.R"))
# source(paste0(root, "core/prob_improv_of_set.R"))
ResolveDependency(root)

dataset <- GetDataReducedAA("SELECT * FROM binary_labeling_activity")
X <- dataset$feature
Y <- dataset$data[, 'type1']
train.X <- X[Y != -1, ]
train.Y <- Y[Y != -1]

# ====================================================================
# Generate recommendation set
print("Generating recommendation set")
gamma.1 <- gamma_1_type1
gamma.0 <- gamma_0_type1
prior.prob <- prior_type1
prior.alpha.1 <- SetPriorReducedAA(gamma.1, NUM_CLASS)
prior.alpha.0 <- SetPriorReducedAA(gamma.0, NUM_CLASS)
# ------------------------------------------------
# Trying a different value
#NUM.RECOM <- 71
#NUM.RECOM <- 20
NUM.RECOM <- 500
# ------------------------------------------------
NB.ITER <- 1000
S <- GenRecomSetOld(train.X, train.Y, prior.alpha.1, prior.alpha.0,
                    prior.prob, NUM.RECOM, NB.ITER, minL, maxL, minR, maxR)


# ====================================================================
# Compute probability improvement
PI.ITER <- 1000
Y <- c(train.Y, rep(0, nrow(S) - 1))
prob.improv.S <- rep(-1, PI.ITER)
rdm.pos.vec <- floor(runif(PI.ITER, 1, nrow(S))) 
print(sprintf("Total Iter = %d", PI.ITER))

for (i in 1:PI.ITER) {
  print(sprintf("iter #%d", i))
  ptm.0 <- proc.time()

  pos <- rdm.pos.vec[i]
  X <- rbind(train.X, S[-pos, ])
  alphas <- BayesianNaiveBayes(X, Y, prior.alpha.1, prior.alpha.0,
                               prior.prob)
  length.left <- ceiling(runif(1, min = minL - 1, max = maxL))
	length.right <- ceiling(runif(1, min = minR - 1, max = maxR))
  best.peptide <- GenOnePeptideMAPOld(length.left, length.right,
                                      alphas$post.alpha.1,
                                      alphas$post.alpha.0, NB.ITER)
  heldout.x <- S[pos, ]
  S[pos, ] <- best.peptide

  ####################################################################
  # test uniqueness of each row of S
  while(nrow(unique(S)) < NUM.RECOM) {
    print("Found a duplicate peptide!")
    length.left <- ceiling(runif(1, min = minL - 1, max = maxL))
	  length.right <- ceiling(runif(1, min = minR - 1, max = maxR))
    best.peptide <- GenOnePeptideMAPOld(length.left, length.right,
                                        alphas$post.alpha.1,
                                        alphas$post.alpha.0, NB.ITER)
    S[pos, ] <- best.peptide
  }
  ####################################################################

  ptm.1 <- proc.time()
  print(sprintf("Until computing PI(S) costs %1.2f", (ptm.1 - ptm.0)[1]))

  prob.improv.S[i] <- ComputeProbImprovOfSet(S, train.X, train.Y,
                                             prior.alpha.1,
                                             prior.alpha.0, 
                                             prior.prob, NB.ITER)

  if (i > 1 && prob.improv.S[i] < prob.improv.S[i - 1]) {
    prob.improv.S[i] <- prob.improv.S[i - 1]
    S[pos, ] <- heldout.x
  }

  print(sprintf("Computing PI(S) costs %1.2f", (proc.time() - ptm.1)[1]))
  print(sprintf("PI(S) = %.3f", prob.improv.S[i]))
  print("==========")
}

