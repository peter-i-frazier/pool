# Local search algorithm
# First run generate_recommendation_MAP_old.R to obtain result$rec
#rm(list=ls())
set.seed(1)

root_path <- ""
datafile <- paste(root_path,"all_TS_data.csv",sep='')

#import libraries
source(paste(root_path,"NB_utility.R",sep=''))
source(paste(root_path,"NB_interface.R",sep=''))
source("prob_improv_of_set.R")

#Specify path/to/classlist file here
classfile <- paste(root_path,"Reduced_AA_Alphabet.csv",sep='')

num_recom <- 71
add_ins <- 1

# In this algorithm we consider type1
gamma_0_type1 <- 50
gamma_1_type1 <- 0.5
prior_type1 <- 1e-4

nF <- 38
nL <- 19
nR <- 19
maxL <- 9
maxR <- 9
minL <- 4
minR <- 4
S.Pos <- nL
# itr <- 500

data_org <- read.csv(datafile, header = T, as.is = T)
AAclass <- read.csv(classfile, header = T, as.is = T)
train_data <- getFeatures(data_org, AAclass, nL, nR)

X_type1 <- train_data[train_data[,'type1'] != -1, 1:(nL+nR)]
Y_type1 <- train_data[train_data[,'type1'] != -1, 'type1']

PI.iter <- 100  # num of iterations running local search
NB.iter <- 1000  # num of iterations running Naive_Bayes
prob.improv.S <- rep(-1, PI.iter)
S <- result$rec
Y <- c(Y_type1, rep(0, dim(S)[1] - 1))
rdn.pos.vec <- floor(runif(PI.iter, 1, dim(S)[1]))  # for randomly selecting a peptide in S to discard from S

print("Ready to run")
print(sprintf("NB.iter = %d", NB.iter))
print(sprintf("Running %d iterations", PI.iter))
# print(rdn.pos.vec)

for (i in 1:PI.iter) {
  print(sprintf("iter #%d", i))
  ptm.0 <- proc.time()

  pos <- rdn.pos.vec[i]
  X <- rbind(X_type1, S[-pos, ])
  alpha <- Dirichlet_Parameter(X, Y, AAclass, gamma_0_type1, gamma_1_type1)
  best.peptide <- old_opt_gen_peptide_lib(nF, maxL, maxR, minL, minR, 
                                          S.Pos, alpha, AAclass)
  heldout.x <- S[pos, ]
  S[pos, ] <- best.peptide

  ptm.1 <- proc.time()
  print(sprintf("Until computing PI(S) costs %1.2f", (ptm.1 - ptm.0)[1]))

  prob.improv.S[i] <- ComputeProbImprovOfSet(S, X_type1, Y_type1, AAclass,
                                             S.Pos, maxL, maxR, 
                                             gamma_0_type1, gamma_1_type1,
                                             prior_type1, NB.iter)
  if (i > 1 && prob.improv.S[i] < prob.improv.S[i - 1]) {
    prob.improv.S[i] <- prob.improv.S[i - 1]
    S[pos, ] <- heldout.x
  }

  print(sprintf("Computing PI(S) costs %1.2f", (proc.time() - ptm.1)[1]))
  print(sprintf("PI(S) = %.3f", prob.improv.S[i]))
  print("==========")
}


