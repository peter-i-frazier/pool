BestProbImprovement.R
	Given training data, computes the probability that at least one peptide
	in a given set of peptides to test is a hit. 

NB.R
	A script using all of the functions to do Naive Bayes training.  Produces output 
	mytest_prob_of_eachAA.csv and prob_CV_reducedAA.csv

NB_predict.R
	Given theta and a prior, predicts the outcome value of each peptide in a set of peptides.

getFeatures.R
	Get the feature vectors (i.e., the reduced amino acid alphabet classes) of all peptides in a set. 

getTheta.R
	A supporting function that computes an estimate of the likelihood parameters
	(theta) of the Naive Bayesian model, assuming the likelihood has (posterior)
	Dirichlet distribution. This estimate is the mean of the Dirichlet
	distribution.

getTheta_MC.R
	A supporting function that simulates the likelihood parameters (theta) of
	the Naive Bayesian model from the posterior, assuming the likelihood has
	(posterior) Dirichlet distribution. 

logistic.R
	Uses logistic regression to do prediction of whether a peptide will be
	a substrate or not.

test/
This directory will contain a bunch of testing files. 
