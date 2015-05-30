# Author: Jialei Wang
# Created: 05.29.2015
# Sample script for generating MAP recommended peptides.

rm(list = ls())
# This section is required in any script, which specifies path of source folder
# and load all modules.
root <- "/Users/jialeiwang/peptide-catalysis/src/"
source(paste0(root, "core/dependency.R"))
ResolveDependency(root)

# Obtaining data, where data is stored in a mysql database. Write mysql query to
# retrieve the data you need. There is a nice tutorial about mysql command here:
# http://www.w3schools.com/sql/default.asp
dataset <- GetDataReducedAA("SELECT * FROM binary_labeling_activity")

# Set prior
alpha.1.prior <- SetPriorReducedAA(0.1, NUM_CLASS)
alpha.0.prior <- SetPriorReducedAA(100, NUM_CLASS)
# Train the naive bayes model
trained.params <- BayesianNaiveBayes(dataset$feature, dataset$data[, 'sfp'],
                                     alpha.1.prior, alpha.0.prior, 1e-4)
# Calculate P(y = 1 | x)
test.peptides.feature <- dataset$feature  # this is just for demonstration
thetas.1 <- SampleThetas(trained.params$post.alpha.1, 1000)
thetas.0 <- SampleThetas(trained.params$post.alpha.0, 1000)
prob.list <- CalculateProb(test.peptides.feature, thetas.1, thetas.0, trained.params$p1)
prob.mean <- prob.list$mean
prob.sd <- prob.list$sd

# Generate peptide recommendation using new method. Here I only show generating one peptide,
# and the recommended peptide is for type 1, i.e., sfp label, AcpS not label, and unlabel.
# Generating the set S is similar to this part of code.
sfp.alpha.1.prior <- SetPriorReducedAA(0.1, NUM_CLASS)
sfp.alpha.0.prior <- SetPriorReducedAA(100., NUM_CLASS)
sfp.p1 <- 1e-4
AcpS.alpha.1.prior <- SetPriorReducedAA(1., NUM_CLASS)
AcpS.alpha.0.prior <- SetPriorReducedAA(100., NUM_CLASS)
AcpS.p1 <- 1e-5
PfAcpH.alpha.1.prior <- SetPriorReducedAA(1., NUM_CLASS)
PfAcpH.alpha.0.prior <- SetPriorReducedAA(50., NUM_CLASS)
PfAcpH.p1 <- 0.5
sfp.trained.params <- BayesianNaiveBayes(dataset$feature, dataset$data[, 'sfp'],
                                         sfp.alpha.1.prior, sfp.alpha.0.prior,
                                         sfp.p1)
AcpS.trained.params <- BayesianNaiveBayes(dataset$feature, dataset$data[, 'AcpS'],
                                          AcpS.alpha.1.prior, AcpS.alpha.0.prior,
                                          AcpS.p1)
PfAcpH.trained.params <- BayesianNaiveBayes(dataset$feature, dataset$data[, 'PfAcpH'],
                                            PfAcpH.alpha.1.prior, PfAcpH.alpha.0.prior,
                                            PfAcpH.p1)
peptide.feature <- GenOnePeptideMAPNew(5, 5, sfp.trained.params$post.alpha.1,
                                       sfp.trained.params$post.alpha.0,
                                       AcpS.trained.params$post.alpha.1,
                                       AcpS.trained.params$post.alpha.0,
                                       PfAcpH.trained.params$post.alpha.1,
                                       PfAcpH.trained.params$post.alpha.0,
                                       1000)

