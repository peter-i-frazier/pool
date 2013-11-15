#This function takes a matrix of peptides as input, each line is a vector with class on each position,
#generate the amino-acid on each position uniformly, and write the generated peptides into a .txt file 
#and a .csv file.

writePep <- function(X, Y, newPep, classlist, Gamma_0, Gamma_1) 
{
	#Compute the probability of being hits of all newPep	
	alpha <- Dirichlet_Parameter(X, Y, classlist, Gamma_0, Gamma_1)
	itr = 1000
	nF <- dim(X)[2]
	nPep <- dim(newPep)[1]
	maxL <- 0
	maxR <- 0
	for (i in 1:(nF/2)) {
		if(any(newPep[,i]!=-1)) {
			maxL = nF/2 - i + 1
			break
		}
	}
	for (i in nF:(nF/2+1) ) {
		if(any(newPep[,i]!=-1)) {
			maxR = i - nF/2
			break
		}
	}
	prob <- 0
	for (i in 1:itr) {
		theta <- getTheta_MC(alpha = alpha, classlist = classlist, Gamma_0 = Gamma_0, Gamma_1 = Gamma_1)
		prob <- NB_predict(newPep, theta, maxL = maxL, maxR = maxR) + prob
	}
	prob <- prob/ itr
	
	#Generate the amino-acid on each position 
	
	
	pepTable <- c()
	for(i in 1:nPep) {
		nTerm <- c()
		for(j in 1:(nF/2)) {
			if( newPep[i,j] != -1 ) {
				AApool <- colnames(classlist)[which(classlist == newPep[i,j])]
				AA <- AApool[sample.int(length(AApool), size = 1)]
				nTerm <- paste(nTerm, AA, sep = "")
			}
		}
		cTerm <- c()
		for(j in (nF/2+1):nF) {
			if( newPep[i,j] != -1 ) {
				AApool <- colnames(classlist)[which(classlist == newPep[i,j])]
				AA <- AApool[sample.int(length(AApool), size = 1)]
				cTerm <- paste(cTerm, AA, sep = "")
			}
		}
		pepTable <- rbind(pepTable, c(nTerm, cTerm))
	}
	pepTable <- cbind(pepTable, prob)
	colnames(pepTable) <- c('nTerm', 'cTerm', 'prob')
	#Sort peptides according to their probabilities of being hits
	pepTable <- pepTable[order(prob, decreasing = TRUE),]
	
	write.csv(pepTable, 'recommendation_AA.csv')
	pepList <- c()
	for( i in 1:nPep ) {
		pep <- paste(pepTable[i,'nTerm'], 'S', pepTable[i, 'cTerm'], sep = "")
		pepList <- rbind(pepList, pep)
	}
	write(pepList, 'recom_list.txt')	
}