#given theta.0  & theta.1
ratio <- theta.1 / theta.0
K <- dim(ratio)[2]
best.class <- rep(0, K)
best.class.ratio <- rep(0,K)
for (i in 1:K) {
	best.class.ratio[i] <- max(ratio[,i])
	best.class[i] <- which(ratio[,i]==max(ratio[,i])
}
# randomly choose length of peptide
L <- ceiling(runif(1)*(maxL - 1)) + 1
R <- ceiling(runif(1)*(maxR - 1)) + 1
best.peptide <- best.class[(19-L+1):(19+R)]