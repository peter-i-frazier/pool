rm(list=ls())
source('gibbs_util.R')
source('../Naive_Bayes/getFeatures.R')
#=================================================================================
#Specify Paths and working directory
dataFile <- '../../data/binaryData_v2.csv'
classFile <- '../../data/Reduced_AA_Alphabet.csv'

#get data 
data.org <- data.frame(read.csv(dataFile, header = T, as.is = T, sep = ","))
AAclass <- read.csv(classFile, header=T, as.is = T, sep=",")

#=================================================================================
#Set parameters
nL <- 19
nR <- 19
trainData <- getFeatures(data.org,AAclass,nL,nR)
nVal <- max(AAclass[1,])
itr <- 1000
Nrec <- 1e6
#=================================================================================
## For AcpH
outcome <- 'AcpH'
X.train <- as.matrix(trainData[,c(1:(nL+nR))])
Y.train <- trainData[,outcome]
Npos <- sum(Y.train)
Nneg <- length(Y.train) - Npos

#calculate mean and std of theta
storage.1 <- matrix(NA, nrow=(nL+nR)*nVal, ncol=itr+2)
storage.0 <- matrix(NA, nrow=(nL+nR)*nVal, ncol=itr+2)
for (n in 1:itr) {
	theta.comb <- sampleTheta(X.train, Y.train, nVal, Npos, Nneg)
	theta.1 <- theta.comb$theta_1
	theta.0 <- theta.comb$theta_0
	for (i in 1:nVal) {
		for (j in 1:(nL+nR)) {
			storage.1[(i-1)*(nL+nR)+j, n] <- theta.1[i,j]
			storage.0[(i-1)*(nL+nR)+j, n] <- theta.0[i,j]
		}
	}
}
for (i in 1:dim(storage.1)[1]) {
	storage.1[i,itr+1] <- mean(storage.1[i, 1:itr])
	storage.1[i,itr+2] <- sd(storage.1[i, 1:itr])
	storage.0[i,itr+1] <- mean(storage.0[i, 1:itr])
	storage.0[i,itr+2] <- sd(storage.0[i, 1:itr])
}
theta.1.sd <- matrix(NA, nrow=nVal, ncol=(nL+nR))
theta.0.sd <- matrix(NA, nrow=nVal, ncol=(nL+nR))
theta.1 <- matrix(NA, nrow=nVal, ncol=(nL+nR))
theta.0 <- matrix(NA, nrow=nVal, ncol=(nL+nR))
for (i in 1:dim(storage.1)[1]) {
	theta.1[ceiling(i/(nL+nR)), i-(nL+nR)*floor((i-1)/(nL+nR))] <- storage.1[i, itr+1]
	theta.1.sd[ceiling(i/(nL+nR)), i-(nL+nR)*floor((i-1)/(nL+nR))] <- storage.1[i, itr+2]
	theta.0[ceiling(i/(nL+nR)), i-(nL+nR)*floor((i-1)/(nL+nR))] <- storage.0[i, itr+1]
	theta.0.sd[ceiling(i/(nL+nR)), i-(nL+nR)*floor((i-1)/(nL+nR))] <- storage.0[i, itr+2]
}
write.csv(theta.1, 'theta_1.csv')
write.csv(theta.0, 'theta_0.csv')
write.csv(theta.1.sd, 'theta_1sd.csv')
write.csv(theta.0.sd, 'theta_0sd.csv')
write.csv(theta.1/theta.0, 'theta_1vstheta_0.csv')

print ('simulate recommendations')
#simulate recommendations
random.peptides <- matrix(NA, nrow=Nrec, ncol=(nL+nR))
random.nL <- ceiling(runif(Nrec) * 3) + 1    # 2-4 uniform
random.nR <- ceiling(runif(Nrec) * 7) + 2    # 3-9 uniform
for (n in 1:Nrec) {
	pp <- runif(random.nL[n]+random.nR[n])
	theta <- theta.1[,(nL-random.nL[n]+1):(nL+random.nR[n])]
	one.peptide <- rep(1,random.nL[n]+random.nR[n])
	for (j in 1:dim(theta)[1]) {
		one.peptide <- one.peptide + ((pp-theta[j,])>0)
		pp <- pp-theta[j,] 
	}
	random.peptides[n,(nL-random.nL[n]+1):(nL+random.nR[n])] <- one.peptide
}
random.peptides <- unique(random.peptides)
prob.random.peptides <- getProb(random.peptides, theta.1, theta.0)
table.forSort <- cbind(prob.random.peptides, random.peptides)
sorted.table <- table.forSort[order(table.forSort[,1],decreasing=T),]
recommend.list <- sorted.table[1:121,-1]
recommend.prob <- sorted.table[1:121,1]

#transfer them into strings & write to file
AA <- c('DE','NQ', 'FWY', 'HKR', 'AILMV', 'GP', 'ST', 'C')
recAAs <- c()
recAA.nterm <- c()
recAA.cterm <- c()
for (i in 1:dim(recommend.list)[1]) {
	one.nterm <- c()
	for (j in 1:nL) {
		if (!is.na(recommend.list[i,j])) {
			AA.group <- unlist(strsplit(AA[recommend.list[i,j]], split=''))
			which.AA <- AA.group[ceiling(length(AA.group) * runif(1))]
			one.nterm <- c(one.nterm, which.AA)
		}
	}
	one.cterm <- c()
	for (j in (nL+1):(nL+nR)) {
		if (!is.na(recommend.list[i,j])) {
			AA.group <- unlist(strsplit(AA[recommend.list[i,j]], split=''))
			which.AA <- AA.group[ceiling(length(AA.group) * runif(1))]
			one.cterm <- c(one.cterm, which.AA)
		}
	}
	recAA.nterm <- c(recAA.nterm, paste(one.nterm, collapse=''))
	recAA.cterm <- c(recAA.cterm, paste(one.cterm, collapse=''))
}
recAA.table <- cbind(recAA.nterm, recAA.cterm)
colnames(recAA.table) <- c('nterm', 'cterm')
write.csv(recAA.table, 'recAA.csv', row.names=F)
# test by looking at recAA.prob and recommend.prob
recAA.prob <- getProb(as.matrix(getFeatures(data.frame(read.csv('recAA.csv', header = T, as.is = T, sep = ",")),AAclass,nL,nR)[,c(1:(nL+nR))]), theta.1, theta.0)

print ('mutate')
# Mutate existing peptides, restrict length to 20
AAlib <- c(colnames(AAclass), 'S')
data_hit <- unique(data.org[data.org$AcpH == 1, ])
nTrain <- dim(data_hit)[1]
recMutates <- c()
recMutates.nterm <- c()
recMutates.cterm <- c()
for (n in 1:1e4) {
	#choose which peptide to mutate
	which.peptide <- ceiling(nTrain * runif(1))
	Nterm <- unlist(strsplit(data_hit[which.peptide,'nterm'], split=''))
	Cterm <- unlist(strsplit(data_hit[which.peptide,'cterm'], split=''))
	L.Nterm <- length(Nterm)
	L.Cterm <- length(Cterm)
	#check if this peptide's length is greater than 20
	if ((L.Nterm+L.Cterm)>19) {
		no.Nterm <- ceiling(7 * runif(1)) + 5
		no.Cterm <- 19 - no.Nterm
		Nterm <- Nterm[(L.Nterm-no.Nterm+1):L.Nterm]
		Cterm <- Cterm[1:no.Cterm]
	}
	L.Nterm <- length(Nterm)
	L.Cterm <- length(Cterm)
	# choose no. of positions to mutate
	no.pos <- ceiling(runif(1)*4)
	# choose which positions to mutate
	which.pos <- ceiling((L.Nterm+L.Cterm) * runif(no.pos))
	#choose AAs to mutate on each position
	for (i in 1:no.pos) {
		if (which.pos[i]<=L.Nterm) {
			while (1) {
				substitute <- AAlib[ceiling(runif(1)*length(AAlib))]
				if (substitute != Nterm[which.pos[i]]) {
					Nterm[which.pos[i]] <- substitute
					break
				}
			}
		} else {
			while (1) {
				substitute <- AAlib[ceiling(runif(1)*length(AAlib))]
				if (substitute != Cterm[which.pos[i]-L.Nterm]) {
					Cterm[which.pos[i]-L.Nterm] <- substitute
					break
				}
			}
		}
	}
	recMutates.nterm <- c(recMutates.nterm, paste(Nterm, collapse=''))
	recMutates.cterm <- c(recMutates.cterm, paste(Cterm, collapse=''))
}
recMutates.table <- cbind(recMutates.nterm, recMutates.cterm)
colnames(recMutates.table) <- c('nterm', 'cterm')
write.csv(recMutates.table, 'recMutates.csv', row.names=F)
# test and sort most probable mutations
recMutates.feature <- as.matrix(getFeatures(data.frame(read.csv('recMutates.csv', header = T, as.is = T, sep = ",")),AAclass,nL,nR)[,c(1:(nL+nR))])
recMutates.table <- recMutates.table[!duplicated(recMutates.feature),]
recMutates.feature <- recMutates.feature[!duplicated(recMutates.feature),]

#test if print 0, it's correct
test.case <- as.matrix(getFeatures(recMutates.table,AAclass,nL,nR))
ttt <- abs(recMutates.feature-test.case)
count <- 0
for (ii in 1:dim(ttt)[1]) {
	for (jj in 1:dim(ttt)[2]) {
		if (!is.na(ttt[ii,jj])) {
			count <- count + ttt[ii,jj]
		}
	}
}
print ('test')
print (count)
###test end

recMutates.prob <- getProb(recMutates.feature, theta.1, theta.0)
sorted.recMutates.table <- recMutates.table[order(recMutates.prob,decreasing=T),]
sorted.recMutates.table <- unique(sorted.recMutates.table)[1:121,]
colnames(sorted.recMutates.table) <- c('nterm', 'cterm')

print ('combine')
# combine recAAs & recMutates by alternating them
REC <- c()
REC.test <- matrix(NA, nrow=121*2, ncol=2)
colnames(REC.test) <- c('nterm', 'cterm')
for (n in 1:121) {
	REC <- c(REC, paste(c(recAA.table[n,1], 'S', recAA.table[n,2]), collapse=''))
	REC <- c(REC, paste(c(sorted.recMutates.table[n,1], 'S', sorted.recMutates.table[n,2]), collapse=''))
	REC.test[n*2-1,] <- recAA.table[n,]
	REC.test[n*2,] <- sorted.recMutates.table[n,]
}
write.csv(REC.test, 'RECtest.csv', row.names=F)
REC.prob <- getProb(as.matrix(getFeatures(data.frame(read.csv('RECtest.csv', header = T, as.is = T, sep = ",")),AAclass,nL,nR)[,c(1:(nL+nR))]), theta.1, theta.0)
# write REC to file
fileRec <- file('recommend.txt')
writeLines(REC, fileRec)
close(fileRec)
