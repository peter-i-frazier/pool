#KG_search: use the knowledge-gradient method to find a set of peptides with 
#a large expected improvement

#nPep: Number of Peptides in the pool
nPep <- 100
#maxIter: Maximum iteration to run the Knowledge-Gradient algorithm
maxIter <- 50
#K: Measure budget for each iteration
K <- 40
#Initialization: Generate nPep peptides randomly
org_peptides <- NB_peptideSim_par(outcome_name = outcome, classlist = AAclass, nPep = 100)
org_peptides <- getFeatures(org_peptides[,c('nterm','cterm')],AAclass,nL,nR)
kg_peptides <- org_peptides
org_proImpr <- ProbImprovement_par(train_AcpH, org_peptides, AAclass, 5000)
org_expImpr <- ExpImprovement_par(train_AcpH, org_peptides, AAclass, 5000)
#Intermediate states storage
Selected <- rep(0,maxIter)
Pos <- rep(0, maxIter)
orgAAclass <- rep(NA, maxIter)
newAAclass <- rep(NA, maxIter)
probI <- rep(NA, maxIter)
expI <- rep(NA, maxIter)

nVal <- length(unique(as.numeric(AAclass)))
for(i in 1:2) {
	#Choose a peptide in the pool and a position in that peptide randomly
	sl_num <- ceiling(nPep*runif(1))
	sl_peptide <- kg_peptides[sl_num,]
	pos <- sample(which(!is.na(sl_peptide)), size = 1)
	#Estimate the measurement precision for each alternative. Currently, we 
	#basically measure each alternative for 10 times and use the inverse of the 
	#sample variance as the measurement precision for each alternative
	m_prec <- rep(0, nVal)
	for(x in 1:nVal) {
		exp_reduc <- rep(0,10)
		sl_peptide[pos] <- x
	    for(j in 1:10){
           	exp_reduc[j] <- ExpImprovement_par(train_AcpH, rbind(kg_peptides[-sl_num,], sl_peptide), AAclass, nRep = 1800)$improve
		}
		m_prec[x] <- 1/var(exp_reduc)
	}
    #Initialize the prior belief of expected improvement resulting from different amino acid class at the selected position
	#kg_theta: The vector of the estimate of the mean of expected improvement resulting from each alternative
	#kg_beta: The vector of measurement precision of each alternative
	#cg_var: The vector of change in the standard deviation of our belief about kg_theta_{n+1} given kg_theta_{n} 
	#kg_theta, kg_beta and cg_var are updated after each run from 1 to N
	kg_theta <- rep(0, nVal)
	kg_beta <- rep(0, nVal)
	cg_error <- rep(Inf, nVal)
	kg_zeta <- rep(0, nVal)
	f <- rep(0, nVal)
	V <- rep(0, nVal)
	obs <- rep(0, K)
	for(n in 1:K) {
		#Compute the knowledge gradient factor for each alternative
		for(x in 1:nVal) {
		    cg_error[x] <- 1/kg_beta[x] - 1/(kg_beta[x]+m_prec[x])	
			kg_zeta[x] <- -abs((kg_theta[x] - max(kg_theta[-x]))/sqrt(cg_error[x]))
			f[x] <- kg_zeta[x]*pnorm(kg_zeta[x])+dnorm(kg_zeta[x])
			V[x] <- sqrt(cg_error[x])*f[x]
		}
		nm <- which.max(V)
		sl_peptide[pos] <- nm
		obs[n] <- ExpImprovement_par(train_AcpH, rbind(kg_peptides[-sl_num,], sl_peptide), AAclass, nRep = 200)$improve
		#Update kg_theta, kg_beta, cg_error for the alternative: nm, using the most recent observation: obs[n]
		kg_theta[nm] <- 1/(kg_beta[nm]+m_prec[nm])*(kg_beta[nm]*kg_theta[nm]+m_prec[nm]*obs[n])		
		kg_beta[nm] <- kg_beta[nm] + m_prec[nm]
        if(kg_theta[nm] == max(kg_theta)) {
			expI[i] <- obs[n] }
	}
	#choose the alternative with largest theta
	orgAAclass[i] <- org_peptides[sl_num,][pos]
    sl_peptide[pos] <- which.max(kg_theta)
	kg_peptides[sl_num,] <- sl_peptide
	#Storage intermediate states:
	Selected[i] <- sl_num
	Pos[i] <- pos
	newAAclass[i] <- sl_peptide[pos]
	probI[i] <- ProbImprovement_par(train_AcpH, kg_peptides, AAclass, 2000)
}
new_proImpr <- ProbImprovement_par(train_AcpH, kg_peptides, AAclass, 5000)
new_expImpr <- ExpImprovement_par(train_AcpH, kg_peptides, AAclass, 5000)