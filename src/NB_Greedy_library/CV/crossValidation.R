args <- commandArgs(TRUE)
for(i in 1:length(args))
{
	eval(parse(text = args[i]))
}

nL <- 19
nR <- 19
S.Pos <- 19
maxL <- 19
maxR <- 19
predIter <- 1000

#Specify path/to/'Naive_Bayes_util.R' here:
source(paste(rootpath,"/src/NB_Greedy_library/Naive_Bayes_util.R",sep=''))

#Specify path/to/classlist file here
classfile <- paste(rootpath,"/data/Reduced_AA_Alphabet.csv",sep='')

data_org <- data.frame(read.csv(datafile, header = T, as.is = T, sep = ','))
AAclass <- data.frame(read.csv(classfile, header = T, as.is = T, sep = ','))
train_data <- getFeatures(data_org, AAclass, nL, nR)
X <- train_data[,1:(nL+nR)]
Y <- train_data[,outcome_name]
print(length(Y))
prob <- Naive_Bayes(X[-n,], Y[-n], X[n,], AAclass, S.Pos, maxL, maxR, Gamma_0, Gamma_1, prior, predIter)
print(prob)
print("gamma 0")
print(Gamma_0)
print("gamma 1")
print(Gamma_1)
print(outcome_name)
write(prob, 'prob')
