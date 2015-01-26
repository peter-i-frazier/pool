root_path = "/fs/home/jw865/peptide-catalysis"
datafile <- paste(root_path,"/data/whole_experiment_data.csv",sep='')

gamma_0_sfp <- 100
gamma_1_sfp <- 0.05
prior_sfp <- 1e-5

gamma_0_AcpS <- 100
gamma_1_AcpS <- 0.1
prior_AcpS <- 0.01

gamma_0_PfAcpH <- 50
gamma_1_PfAcpH <- 0.5
prior_PfAcpH <- 0.1

# gamma_0_type1 <- 50
# gamma_1_type1 <- 0.5
# prior_type1 <- 1e-4
# 
# gamma_0_type2 <- 1000
# gamma_1_type2 <- 0.5
# prior_type2 <- 1e-4

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

X_PfAcpH <- train_data[train_data[,'PfAcpH'] != -1,1:(nL+nR)]
Y_PfAcpH <- train_data[train_data[,'PfAcpH'] != -1, 'PfAcpH']
X_sfp <- train_data[train_data[,'sfp'] != -1, 1:(nL+nR)]
Y_sfp <- train_data[train_data[,'sfp'] != -1, 'sfp']
X_AcpS <- train_data[train_data[,'AcpS'] != -1, 1:(nL+nR)]
Y_AcpS <- train_data[train_data[,'AcpS'] != -1, 'AcpS']
X_type1 <- train_data[train_data[,'type1'] != -1, 1:(nL+nR)]
Y_type1 <- train_data[train_data[,'type1'] != -1, 'type1']
X_type2 <- train_data[train_data[,'type2'] != -1, 1:(nL+nR)]
Y_type2 <- train_data[train_data[,'type2'] != -1, 'type2']


table <- c()
theta <- NB_theta_matrices(X_sfp, Y_sfp, AAclass, S.Pos, maxL, maxR, gamma_0_sfp, gamma_1_sfp, prior_sfp)
write.csv(theta$theta_1 / theta$theta_0, "NB_illustration.csv")
