#Need to Run 'NB_Simulation.R' first to import data!
#===========================================================================
#Test Case 2
#Generating true theta
theta_true <- getTheta_MC(alpha = alpha, nVal = nVal)

par(mfrow = c(2,3))
prior.positive <- c(10**-(4:1),1/2)
for (i in 1:length(prior.positive)) {
    #Simulate training data
    nTrain <- 10**4
    train_1 <- NB_peptideSim(theta = theta_true, outcome_name = outcome, nVal = nVal, nPep = nTrain/2, isHit = 1)
    train_0 <- NB_peptideSim(theta = theta_true, outcome_name = outcome, nVal = nVal, nPep = nTrain/2, isHit = 0)
    sTrain <- rbind(train_1, train_0)
	sTrainData <- getFeatures(sTrain, class ,nL, nR)
    #Generate test data using simulated training data
    nTestPep <- 100
    sTestData <- NB_peptideSim(nVal = nVal, nPep = nTestPep)
    sTestFeat <- getFeatures(sTestData[,c('nterm','cterm')], class, nL, nR)
    #Compute probability each test peptide is a hit using "true" theta
    prob_trueTheta <- NB_predict(sTestFeat, theta_true, prior.positive = prior.positive[i])
	#Compute probability each test peptide is a hit using simulation
	s_predict_mat <- c()
	s_alpha <- Dirichlet_Parameter(sTrainData,nVal) 
	nRep <- 100
    for(j in 1:nRep) {
    sTheta <- getTheta_MC(alpha = s_alpha, nVal = nVal)
    s_predict_mat <- rbind(s_predict_mat,
        NB_predict(sTestFeat, sTheta, prior.positive = prior.positive[i]))
    }
    prob_simuTheta <- colMeans(s_predict_mat)
	#Plot prob_trueTheta vs. prob_simuTheta
    plot(prob_trueTheta, prob_simuTheta)
}

