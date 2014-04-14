#Specify path/to/'Naive_Bayes_util.R' here:
source('/fs/home/py75/Documents/peptide/src/NB_Greedy_library/Naive_Bayes_util.R')

#Specify path/to/training data file here
datafile <- 'fs/home/py75/Documents/peptide/data/Data_10.csv'

#Specify path/to/classlist file here
classfile <- 'fs/home/py75/Documents/peptide/data/Reduced_AA_Alphabet.csv'

args <- commandArgs(TRUE)
print(length(args))

for(i in 1:length(args)) {
	eval(parse(text = args[i]))
	print(args[i])
}
data_org <- data.frame(read.csv(datafile, header = T, as.is = T, sep = ','))
AAclass <- data.frame(read.csv(classfile, header = T, as.is = T, sep = ','))
train_data <- getFeatures(data_org, AAclass, nL, nR)
X <- train_data[,1:(nL+nR)]
Y <- train_data[,outcome_name]
output <- Naive_Bayes(X[-n,], Y[-n], X[n,], AAclass, S.Pos, maxL, maxR, Gamma_0, Gamma_1, prior.positive, predIter)
write(output, 'output')
