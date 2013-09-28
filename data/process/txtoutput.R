rm(list=ls())
source('../../src/Naive_Bayes/Naive_Bayes_util.R')
source('../../src/Opt_Search/mutateREC.R')
#=================================================================================
#Specify Paths and working directory
dataFile <- 'Data.csv'
classFile <- '../Reduced_AA_Alphabet.csv'

#get data 
data.org <- data.frame(read.csv(dataFile, header = T, as.is = T, sep = ","))
# AA -> class  & class -> AA
classlist <- read.csv(classFile, header=T, as.is = T, sep=",")

data <- getFeatures(data.org, classlist, 19, 19)
data <- data[,-39]
mutate.result <- mutateREC(data, classlist, 150)
mutate <- cbind(mutate.result$rec, mutate.result$prob)
write.csv(mutate, 'mutation.csv', row.names=F)


recom_1 <- as.matrix(read.csv('recom_1_sort.csv'))
recom_2 <- as.matrix(read.csv('recom_2_sort.csv'))
recom_3 <- as.matrix(read.csv('recom_3_sort.csv'))
recom_4 <- as.matrix(read.csv('recom_4_sort.csv'))

# recom_1 <- recom_1[order(recom_1[,39],decreasing=T),]
# recom_2 <- recom_2[order(recom_2[,39],decreasing=T),]
# recom_3 <- recom_3[order(recom_3[,39],decreasing=T),]
# recom_4 <- recom_4[order(recom_4[,39],decreasing=T),]

# write.csv(recom_1, 'recom_1_sort.csv', row.names=F)
# write.csv(recom_2, 'recom_2_sort.csv', row.names=F)
# write.csv(recom_3, 'recom_3_sort.csv', row.names=F)
# write.csv(recom_4, 'recom_4_sort.csv', row.names=F)
mix_idx <- (1:10) * 5
mix <- rbind(recom_1[mix_idx,], recom_2[mix_idx,], recom_3[mix_idx,], recom_4[mix_idx,])
recom_1 <- recom_1[-mix_idx,]
recom_2 <- recom_2[-mix_idx,]
recom_3 <- recom_3[-mix_idx,]
recom_4 <- recom_4[-mix_idx,]

AA <- c('DE','NQ', 'FWY', 'HKR', 'AILMV', 'GP', 'ST', 'C')
feature.table <- rbind(recom_1, recom_2, recom_3, recom_4, mutate)
AA.table <- c()
for (i in 1:dim(feature.table)[1]) {
	feature <- feature.table[i,-39]
	nterm <- ''
	for (j in 1:19) {
		if (feature[j] != -1) {
			class <- unlist(strsplit(AA[feature[j]], split=''))
			nterm <- paste0(nterm, class[ceiling(length(class)*runif(1))])
		}
	}
	cterm <- ''
	for (j in 20:38) {
		if (feature[j] != -1) {
			class <- unlist(strsplit(AA[feature[j]], split=''))
			cterm <- paste0(cterm, class[ceiling(length(class)*runif(1))])
		}
	}
	AA.table <- rbind(AA.table, c(nterm, cterm, feature.table[i,39]))
}
colnames(AA.table) <- c('nterm', 'cterm', 'prob')

recom_1_AA <- AA.table[1:140,]
recom_2_AA <- AA.table[141:280,]
recom_3_AA <- AA.table[281:420,]
recom_4_AA <- AA.table[421:560,]
mutate_AA <- AA.table[561:710,]
write.csv(recom_1_AA, 'recom_1_AA.csv', row.names=F)
write.csv(recom_2_AA, 'recom_2_AA.csv', row.names=F)
write.csv(recom_3_AA, 'recom_3_AA.csv', row.names=F)
write.csv(recom_4_AA, 'recom_4_AA.csv', row.names=F)
write.csv(mutate_AA, 'mutate_AA.csv', row.names=F)

# deal with mix data
reference <- unlist(strsplit('zzzzzzzNzASFzzDLGADLDTzELVzzzzzzzzzzzz', split=''))
ref.class <- c(0,0,0,0,0,0,0,2,0,5,7,3,0,0,1,5,6,5,1,5,1,7,0,1,5,5,0,0,0,0,0,0,0,0,0,0,0,0)
mix_AA <- c()
for (i in 1:40) {
	feature <- mix[i,-39]
	compare <- feature - ref.class
	nterm <- ''
	for (j in 1:19) {
		if (compare[j] == 0) {
			nterm <- paste0(nterm, reference[j])
		} else if (feature[j] != -1) {
			class <- unlist(strsplit(AA[feature[j]], split=''))
			nterm <- paste0(nterm, class[ceiling(length(class)*runif(1))])
		}
	}
	cterm <- ''
	for (j in 20:38) {
		if (compare[j] == 0) {
			cterm <- paste0(cterm, reference[j])
		} else if (feature[j] != -1) {
			class <- unlist(strsplit(AA[feature[j]], split=''))
			cterm <- paste0(cterm, class[ceiling(length(class)*runif(1))])
		}
	}
	mix_AA <- rbind(mix_AA, c(nterm, cterm, mix[i,39]))
}
colnames(mix_AA) <- c('nterm', 'cterm', 'prob')
write.csv(mix_AA, 'mix_AA.csv', row.names=F)

REC <- c()
for (i in 1:40) {
	REC <- c(REC, paste0(recom_1_AA[i,'nterm'], 'S', recom_1_AA[i,'cterm']))
	REC <- c(REC, paste0(recom_2_AA[i,'nterm'], 'S', recom_2_AA[i,'cterm']))
	REC <- c(REC, paste0(recom_3_AA[i,'nterm'], 'S', recom_3_AA[i,'cterm']))
	REC <- c(REC, paste0(recom_4_AA[i,'nterm'], 'S', recom_4_AA[i,'cterm']))
	REC <- c(REC, paste0(mutate_AA[i,'nterm'], 'S', mutate_AA[i,'cterm']))
	REC <- c(REC, paste0(mix_AA[i,'nterm'], 'S', mix_AA[i,'cterm']))
}
for (i in 41:140) {
	REC <- c(REC, paste0(recom_1_AA[i,'nterm'], 'S', recom_1_AA[i,'cterm']))
	REC <- c(REC, paste0(recom_2_AA[i,'nterm'], 'S', recom_2_AA[i,'cterm']))
	REC <- c(REC, paste0(mutate_AA[i,'nterm'], 'S', mutate_AA[i,'cterm']))
}
for (i in 141:150) {
	REC <- c(REC, paste0(mutate_AA[i,'nterm'], 'S', mutate_AA[i,'cterm']))
}

# test REC has no lines greater than 20
isbad <- 0
for (i in 1:length(REC)) {
	pep <- unlist(strsplit(REC[i], split=''))
	if (length(pep) > 20) {
		isbad <- 1
		break
	}
}
print (paste('greater than 20', isbad))

# write REC to file
fileRec <- file('recommendation.txt')
writeLines(REC, fileRec)
close(fileRec)