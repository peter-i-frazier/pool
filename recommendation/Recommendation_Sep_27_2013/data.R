mutate <- read.csv('mutate_AA.csv',as.is=T)
data <- read.csv('~/peptide-catalysis/data/newData#2.csv', as.is=T)
data <- cbind(data, rep(0,dim(data)[1]))
for (n in 1:dim(data)[1]) {
    for (i in 1:dim(mutate)[1]) {
        if (data[n,'nterm'] == mutate[i,'nterm'] && data[n,'cterm'] == mutate[i,'cterm']) {
            data[n,5] = 1
            break
        }
    }
}
write.csv(data, 'modified_data.csv')
