rm(list=ls())

source('gibbs_util.R')

feature <- as.matrix(read.csv('feature.csv', header=F, as.is=T))
feature <- feature[-1,]

# leave one out cross validation
prob.hit <- c()
for (N in 1:8) {
	print (N)
	# Initialization
	train.data <- feature[1:8,]
	test.data <- train.data[N,]
	train.data <- train.data[-N,]
	prob.1 <- gibbsSampler(train.data, test.data)
	prob.hit <- c(prob.hit, prob.1)
}

# take in all the data and predict prob not hit
train.data <- feature[1:8,]
test.data <- feature[9:13,]
prob.0 <- gibbsSampler(train.data, test.data)

write.csv(prob.hit,'result/prob_hit_enzyme2.csv')
write.csv(prob.0,'result/prob_nothit_enzyme2.csv')
