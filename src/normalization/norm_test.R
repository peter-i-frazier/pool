rm(list = ls())
source('norm_util.R')

simulate_data <- function(mu, sigma, num_data, j_class) {
    seed <- read.csv('data/simulate_seed.csv', stringsAsFactors=F)
    seed <- seed[1:num_data,]
    theta <- seed[, 'sfp_1']
    simulate_data <- c()
    for (i in 1:dim(seed)[1]) {
        num_rep <- sample(2:4, 1)
        for (j in sample(1:length(j_class), num_rep)) {
            to_bind <- seed[i,]
            to_bind['TS'] <- j_class[j]
            to_bind['sfp_1'] <- (theta[i] + rnorm(1,0,0)) * sigma[j] + mu[j]
            #to_bind['sfp_1'] <- (theta[i]) * sigma[j] + mu[j]
            simulate_data <- rbind(simulate_data, to_bind)
        }
    }
    return (list(theta=theta, data=simulate_data, sum=sum(theta)))
}

group <- 1
if (group == 1) {
    mu <- rep(0, 5)
    sigma <- c(1,4,10,15,1)
    num_data <- 30
    j.class <- 1:5
} else {
    mu <- rep(0, 4)
    sigma <- c(1,4,10,15)
    num_data <- 30
    j.class <- c(1,3,4,5)
}

eps <- 1e-5
simu.list <- simulate_data(mu, sigma, num_data, j.class)
reduced.data <- simu.list$data
unique.seqs <- unique(reduced.data[, 'seq'])
dup.data <- c()
for (seq in unique.seqs) {
    sub.data <- reduced.data[reduced.data[, 'seq'] == seq, c('TS', 'spot', 'seq', 'sfp_1')]
    if (nrow(sub.data) > 1) {
        dup.data <- rbind(dup.data, sub.data)
    }
}
sol <- normalize(reduced.data, dup.data, 'sfp_1', const=simu.list$sum)
print (ifelse(sum(sigma - 1 / sol[1:length(j.class)]) < eps, 'pass', 'fail'))
