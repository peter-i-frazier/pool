# seq_idx <- 1
# root_path <- "/fs/home/jw865/peptide-catalysis"
# datafile <- paste(root_path,"/data/all_TS_data.csv",sep='')

args <- commandArgs(TRUE)
for(i in 1:length(args))
{
	eval(parse(text = args[i]))
}

gamma_0_sfp <- 100
gamma_1_sfp <- 0.01
prior_sfp <- 1e-4

gamma_0_AcpS <- 100
gamma_1_AcpS <- 0.1
prior_AcpS <- 1e-4

gamma_0_PfAcpH <- 1000
gamma_1_PfAcpH <- 10
prior_PfAcpH <- 1e-4

nL <- 19
nR <- 19
S.Pos <- nL
maxL <- 19
maxR <- 19
itr <- 500

#import libraries
source(paste(root_path,"/src/NB_Greedy_library/NB_utility.R",sep=''))
source(paste(root_path,"/src/NB_Greedy_library/NB_interface.R",sep=''))

#Specify path/to/classlist file here
classfile <- paste(root_path,"/data/Reduced_AA_Alphabet.csv",sep='')

data_org <- read.csv(datafile, header = T, as.is = T)
AAclass <- read.csv(classfile, header = T, as.is = T)
train_data <- getFeatures(data_org, AAclass, nL, nR)

predict_X <- train_data[seq_idx, 1:(nL+nR)]
reduced_train_data <- train_data[-seq_idx,]

X_PfAcpH <- reduced_train_data[reduced_train_data[,'PfAcpH'] != -1,1:(nL+nR)]
Y_PfAcpH <- reduced_train_data[reduced_train_data[,'PfAcpH'] != -1, 'PfAcpH']
X_sfp <- reduced_train_data[reduced_train_data[,'sfp'] != -1, 1:(nL+nR)]
Y_sfp <- reduced_train_data[reduced_train_data[,'sfp'] != -1, 'sfp']
X_AcpS <- reduced_train_data[reduced_train_data[,'AcpS'] != -1, 1:(nL+nR)]
Y_AcpS <- reduced_train_data[reduced_train_data[,'AcpS'] != -1, 'AcpS']


# new method
prob_sfp <- Naive_Bayes(X_sfp, Y_sfp, predict_X, AAclass, S.Pos, nL, nR, gamma_0_sfp, gamma_1_sfp, prior_sfp, itr)
prob_AcpS <- Naive_Bayes(X_AcpS, Y_AcpS, predict_X, AAclass, S.Pos, nL, nR, gamma_0_AcpS, gamma_1_AcpS, prior_AcpS, itr)
prob_PfAcpH<- Naive_Bayes(X_PfAcpH, Y_PfAcpH, predict_X, AAclass, S.Pos, nL, nR, gamma_0_PfAcpH, gamma_1_PfAcpH, prior_PfAcpH, itr)

write.csv(c(prob_sfp, prob_AcpS, prob_PfAcpH, train_data[seq_idx, 'org_idx']), paste('prob_', seq_idx, '.csv', sep='')) 
