I <- 30
J <- 12
K <- 8
a <- round(runif(I) * 10,3)
ita <- array(0, dim=c(I,K,J))
for (i in 1:I) {
	for (j in 1:J) {
		for (k in 1:K) {
			ita[i,k,j] <- round(runif(1)*2,3)
		}
	}
}

# write to file
f <- file('test_data_I30.dat','a')
idxI <- rep(0,I)
for (i in 1:I) {
	idxI[i] <- paste(c('I',i), collapse="")
}
writeLines(paste(c('set I := ', idxI, ';'), collapse=' '),f)

idxJ <- rep(0,J)
for (j in 1:J) {
	idxJ[j] <- paste(c('J',j), collapse="")
}
writeLines(paste(c('set J := ', idxJ, ';'), collapse=' '),f)

idxK <- rep(0,K)
for (k in 1:K) {
	idxK[k] <- paste(c('K',k), collapse="")
}
writeLines(paste(c('set K := ', idxK, ';'), collapse=' '),f)
writeLines('\n', f)

str_a <- rep(0,I)
for (i in 1:I) {
	str_a[i] <- paste(c(idxI[i], a[i]), collapse=' ')
}
writeLines(paste(c('param a := ', str_a, ';'), collapse=' '),f)
writeLines('\n', f)

writeLines('param ita :=', f)
for (i in 1:I) {
	writeLines('\n', f)
	head <- paste(c("['", idxI[i], "',*,*]:"), collapse='')
	writeLines(paste(c(head, idxJ, ':='), collapse=' '), f)
	for (k in 1:K) {
		if (i==I && k==K) {
			writeLines(paste(c(idxK[k], '           ', ita[i,k,], ';'), collapse=' '), f)
		} else {
			writeLines(paste(c(idxK[k], '           ', ita[i,k,]), collapse=' '), f)
		}
	}
}
close(f)