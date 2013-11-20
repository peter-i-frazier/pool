load('rec.Rdata')
N <- 100
b <- 12
P.recom <- rep(0,N)
P.naive <- rep(0,N)
P.mutate <- rep(0,N)
for (i in 1:N) {
	P.recom[i] <- prob_shortest_hit(recom$rec[1:i,], recom_prob[1:i], b)
	P.naive[i] <- prob_shortest_hit(naive_rec[1:i,], naive_prob[1:i], b)
	
}
