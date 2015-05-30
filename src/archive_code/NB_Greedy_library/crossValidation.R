args <- commandArgs(TRUE)
for(i in 1:length(args))
{
	eval(parse(text = args[i]))
}

gamma0 <- c(1, 10, 100)
gamma1 <- c(0.1, 0.5, 1)
Gamma_0 <- gamma0[ceiling(gamma_idx/length(gamma1))]
Gamma_1 <- gamma1[gamma_idx - floor(gamma_idx/length(gamma1))*length(gamma1)]
ROOT_PATH = "/fs/home/jw865/peptide-catalysis"

#Specify path/to/'Naive_Bayes_util.R' here:
source(paste(ROOT_PATH,"/src/NB_Greedy_library/Naive_Bayes_util.R",sep=''))
#Specify path/to/training data file here
datafile <- paste(ROOT_PATH,"/data/2014_06_03_orthogonal_labeling_data/Training_Set_Cumulative.csv",sep='')

#Specify path/to/classlist file here
classfile <- paste(ROOT_PATH,"/data/Reduced_AA_Alphabet.csv",sep='')

data_org <- data.frame(read.csv(datafile, header = T, as.is = T, sep = ','))
AAclass <- data.frame(read.csv(classfile, header = T, as.is = T, sep = ','))
train_data <- getFeatures(data_org, AAclass, nL, nR)
X <- train_data[,1:(nL+nR)]
Y <- train_data[,outcome_name]
print(length(Y))
output <- Naive_Bayes(X[-n,], Y[-n], X[n,], AAclass, S.Pos, maxL, maxR, Gamma_0, Gamma_1, prior.positive, predIter)
print(output)
print("gamma 0")
print(Gamma_0)
print("gamma 1")
print(Gamma_1)
write(output, 'output')
