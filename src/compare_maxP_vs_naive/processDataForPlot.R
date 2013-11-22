rm(list=ls())
source('prob_shortest_hit.R')
load('recom_data#2Benchmark.RData')
N <- 100
b <- 12
P.recom2 <- rep(0,N)
P.naive <- rep(0,N)
P.mutate <- rep(0,N)
for (i in 1:N) {
	P.recom2[i] <- prob_shortest_hit(recom2[1:i,], recom2_prob[1:i], b)
	P.naive[i] <- prob_shortest_hit(naive_rec[1:i,], naive_prob[1:i], b)
	P.mutate[i] <- prob_shortest_hit(mutate_rec[1:i,], mutate_prob[1:i], b)
}
# save(list= c('P.recom2', 'P.naive', 'P.mutate'), file='P for plot')
PP <- rbind(P.recom2, P.naive, P.mutate)
write.csv(PP, 'PP2.csv')