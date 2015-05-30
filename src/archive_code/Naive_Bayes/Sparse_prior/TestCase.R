#A test script for the Sparse Prior approach


#Step 1: Simulate the true theta
set.seed(1)
source(paste(srcPath, '/getTheta_MC.R',sep=""))
trainData <- trainData[,c(1:(nL+nR),which(colnames(trainData)==outcome))]
theta_true <- getTheta_MC(train = trainData, classlist = AAclass)

#Step 2:Simulate Training Data
source(paste(srcPath, '/NB_peptideSim_par.R', sep=""))
nTrain <- 10**4
train_1 <- NB_peptideSim_par(theta = theta_true, outcome_name = outcome, classlist = AAclass, nPep = nTrain/2, isHit = 1, minL = 9, maxL = 9, minR = 9, maxR = 9)
rownames(train_1) <- NULL
train_0 <- NB_peptideSim_par(theta = theta_true, outcome_name = outcome, classlist = AAclass, nPep = nTrain/2, isHit = 0, minL = 9, maxL = 9, minR = 9, maxR = 9)
rownames(train_0) <- NULL
sTrain <- rbind(train_1, train_0)
sTrain <- getFeatures(sTrain, AAclass ,nL, nR)
sTrain.X <- sTrain[,c(1:(nL+nR))]
sTrain.Y <- sTrain[,outcome]

#Step 3:Simulate Test Data
nTestPep <- 100
sTest <- NB_peptideSim_par(classlist = AAclass, nPep = nTestPep, minL = 9, maxL = 9, minR = 9, maxR = 9)
sTest <- getFeatures(sTest[,c('nterm','cterm')], AAclass, nL, nR)

#Step 4:Predict probability of being hits of test data using true theta
source(paste(srcPath,'/NB_predict_par.R',sep=""))
prob_trueTheta <- NB_predict_par(sTest, theta_true, prior.positive = 1e-1) 

#Step 5 Predict probability of being hits of test data using training data
prob_trainData <- sparsePrior(sTrain.X, sTrain.Y, sTest, nAA, burnin.step, record.step)