rm(list=ls())

root_path <- "/fs/home/jw865/peptide-catalysis"
datafile <- paste(root_path,"/data/whole_experiment_data.csv",sep='')
#import libraries
source(paste(root_path,"/src/NB_Greedy_library/NB_utility.R",sep=''))
source(paste(root_path,"/src/NB_Greedy_library/NB_interface.R",sep=''))

#Specify path/to/classlist file here
classfile <- paste(root_path,"/data/Reduced_AA_Alphabet.csv",sep='')

# parameters that can change
num_recom <- 50

gamma_0_sfp <- 100
gamma_1_sfp <- 0.5
prior_sfp <- 1e-4

gamma_0_AcpS <- 100
gamma_1_AcpS <- 0.5
prior_AcpS <- 1e-4

gamma_0_PfAcpH <- 10
gamma_1_PfAcpH <- 0.05
prior_PfAcpH <- 0.5

gamma_0_type1 <- 100
gamma_1_type1 <- 0.01
prior_type1 <- 1e-5

gamma_0_type2 <- 10
gamma_1_type2 <- 1
prior_type2 <- 1e-5

add_ins <- 1
nL <- 19
nR <- 19
maxL <- 9
maxR <- 9
minL <- 4
minR <- 4
S.Pos <- nL
itr <- 500
# end parameters

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

# # find type 1 hit new method
# method <- paste('type_1_new_addins_', toString(add_ins), sep = '')
# result <- generate_recommendation_MAP_new(X_sfp, Y_sfp, X_AcpS, Y_AcpS, X_PfAcpH, Y_PfAcpH, AAclass, S.Pos, num_recom, maxL, maxR, minL, minR, gamma_0_sfp, gamma_1_sfp, prior_sfp, gamma_0_AcpS, gamma_1_AcpS, prior_AcpS, gamma_0_PfAcpH, gamma_1_PfAcpH, prior_PfAcpH, add_ins, itr)

# # find type 2 hit new method
# method <- paste('type_2_new_addins_', toString(add_ins), sep = '')
# result <- generate_recommendation_MAP_new(X_AcpS, Y_AcpS, X_sfp, Y_sfp, X_PfAcpH, Y_PfAcpH, AAclass, S.Pos, num_recom, maxL, maxR, minL, minR, gamma_0_AcpS, gamma_1_AcpS, prior_AcpS, gamma_0_sfp, gamma_1_sfp, prior_sfp, gamma_0_PfAcpH, gamma_1_PfAcpH, prior_PfAcpH, add_ins, itr)

# find type 1 hit old method
method <- paste('type_1_old_addins_', toString(add_ins), sep = '')
result <- generate_recommendation_MAP_old(X_type1, Y_type1, AAclass, S.Pos, num_recom, maxL, maxR, minL, minR, gamma_0_type1, gamma_1_type1, prior_type1, add_ins, itr)

# # find type 2 hit old method
# method <- paste('type_2_old_addins_', toString(add_ins), sep = '')
# result <- generate_recommendation_MAP_old(X_type2, Y_type2, AAclass, S.Pos, num_recom, maxL, maxR, minL, minR, gamma_0_type2, gamma_1_type2, prior_type2, add_ins, itr)

writePep(result$rec, S.Pos, result$prob, AAclass, paste(method, '.csv', sep=''), method)
print (method)
print (add_ins)


