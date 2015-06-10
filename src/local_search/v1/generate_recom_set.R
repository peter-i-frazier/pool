# Author: Tom Fei
# Created: 06.06.2015
# Generate a recommendation set
# Use type1 response

rm(list = ls())
set.seed(1)
# This section is required in any script, which specifies path of source folder
# and load all modules.
root <- "~/repos/peptide-catalysis/src/"
source(paste0(root, "core/dependency.R"))
ResolveDependency(root)

dataset <- GetDataReducedAA("SELECT * FROM binary_labeling_activity")
X <- dataset$feature
Y <- dataset$data[, 'type1']
X <- X[Y != -1, ]
Y <- Y[Y != -1]
NUM.RECOM <- 71
NB.ITER <- 1000
num.unique.recom <- 0
recom.set <- c()

gamma.1 <- gamma_1_type1
gamma.0 <- gamma_0_type1
prior.prob <- prior_type1
prior.alpha.1 <- SetPriorReducedAA(gamma.1, NUM_CLASS)
prior.alpha.0 <- SetPriorReducedAA(gamma.0, NUM_CLASS)

# length.left <- 6
# length.right <- 6

print("Start recommending")
#for (i in 1:num.recom) {
while (num.unique.recom < NUM.RECOM) { 
  print(sprintf("iter %d", num.unique.recom))

  alphas <- BayesianNaiveBayes(X, Y, prior.alpha.1, prior.alpha.0,
                               prior.prob)
  
  #print("Passed NB")

  length.left <- ceiling(runif(1, min = minL - 1, max = maxL))
	length.right <- ceiling(runif(1, min = minR - 1, max = maxR))
  #length.left <- sample(minL:maxL, 1)
  #length.right <- sample(minR:maxR, 1)
  new.rec.peptide <- GenOnePeptideMAPOld(length.left, length.right,
                                         alphas$post.alpha.1,
                                         alphas$post.alpha.0, NB.ITER)

  #print(sprintf("row = %d", nrow(unique(rbind(recom.set, new.rec.peptide)))))
  if (nrow(unique(rbind(recom.set, new.rec.peptide))) ==
      num.unique.recom)  # already contained this peptide
    next
  num.unique.recom <- num.unique.recom + 1
  X <- rbind(X, new.rec.peptide)
  Y <- c(Y, 0)
  recom.set <- rbind(recom.set, new.rec.peptide)
}

