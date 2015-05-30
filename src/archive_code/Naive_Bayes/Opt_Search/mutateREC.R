#================================================================================
#Function: mutateREC
#
#--------------------------------------------------------------------------------
#Description:
#    Mutate existing data and give the first 150 peptides with highest prob being hit
#
#--------------------------------------------------------------------------------
#Input arguments:
#	 data
#			 data.frame with first 38 columns as features & last column response
#	 classlist
#             A matrix specifying the mapping between each amino-acids and the 
#             class it belongs to.
#    Nrec
# 			 No. of recommendations to be generated
# 	 itr
#			 Iteration needed for calculation prob being hit
#	 Nlib
# 			 No. of peptides considered for searching optimal
#	
#
#--------------------------------------------------------------------------------
#Return objects:
#    a list
#            result$rec:	matrix with each row representing a recommended peptide  
#--------------------------------------------------------------------------------
mutateREC <- function(data, classlist, Nrec, itr=500, Nlib = 1e5) {
	nAA <- length(unique(as.numeric(classlist)))
	Y <- data[,dim(data)[2]]
	class(Y)
	X <- as.matrix(data[,-dim(data)[2]])
	Ntrain <- dim(X)[1]

	rec.feature <- matrix(-1, nrow=Nlib, ncol = 38)
	which.peptide <- ceiling(Ntrain * runif(Nlib))   # choose which peptide to mutate
	no.pos <- ceiling(runif(Nlib)*4)  # choose no. of positions to mutate
	for (n in 1:Nlib) {
		print (n)
		to.mutate <- X[which.peptide[n],]
		# locate start and end point
		for (i in 1:38) {
			if (to.mutate[i] != -1) {
				start <- i
				break
			}
		}
		for (i in start:38) {
			if (to.mutate[i] == -1) {
				end <- i-1
				break
			}
		}
		# locate mutating position
		pos <- ceiling(runif(no.pos[n]) * (end-start+1)) + start - 1
		for (k in 1:length(pos)) {
			to.mutate[pos[k]] <- ceiling(runif(1) * nAA)
		}
		rec.feature[n,] <- to.mutate
	}
	# calculate prob
	theta.1 <- 0
	theta.0 <- 0
	for (i in 1:itr) {
			print (i)
			theta <- getTheta_MC(X, Y, classlist=classlist, Gamma = 0.3)
			theta.1 <- theta.1 + theta$theta_1
			theta.0 <- theta.0 + theta$theta_0
		}
	theta.1 <- theta.1 / itr
	theta.0 <- theta.0 / itr
	theta <- list('theta_0'=theta.0, 'theta_1'=theta.1)
	prob <- NB_predict(rec.feature, theta)
	# delete duplicate and choose first Nrec peptides
	sorted_rec.feature <- rec.feature[order(prob, decreasing=T),]
	prob <- prob[order(prob, decreasing=T)]
	final.rec <- matrix(-1, nrow=Nrec, ncol=38)
	final.prob <- c()
	count <- 1
	for (i in 1:dim(sorted_rec.feature)[1]) {
		print (i)
		print (count)
		if (count==1) {
			duplicate <- 0
			for (j in 1:dim(X)[1]) {
				if (prod(X[j,] == sorted_rec.feature[i,])==1) {
					duplicate <- 1
					break
				}
			}
			if (duplicate == 0) {
				final.rec[count,] <- sorted_rec.feature[i,]
				final.prob <- c(final.prob, prob[i])
				count <- count +1
			} 
		} else if (prod(sorted_rec.feature[i,]==sorted_rec.feature[i-1,])==0) {
			duplicate <- 0
			for (j in 1:dim(X)[1]) {
				if (prod(X[j,] == sorted_rec.feature[i,])==1) {
					duplicate <- 1
					break
				}
			}
			if (duplicate == 0) {
				final.rec[count,] <- sorted_rec.feature[i,]
				final.prob <- c(final.prob, prob[i])
				count <- count +1
			}
		}
		if (count > Nrec) {
			break
		}
	}
	result <- list('rec'=final.rec, 'prob'=final.prob)
	return (result)
}