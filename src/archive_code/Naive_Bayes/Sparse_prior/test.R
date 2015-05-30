# # burn in step
# for (t in 1:3000) {
# 	print (t)
# 	mu.1 <- getMu(X.1, factor, Z.1, P.1)
# 	mu.0 <- getMu(X.0, 1, Z.0, P.0)
# 	Z.1 <- getZ(X.1, P.1, mu.1, ztable)
# 	Z.0 <- getZ(X.0, P.0, mu.0, ztable)
# 	P.1 <- getP(Z.1)
# 	P.0 <- getP(Z.0)
# }
#test stability 1
mumu.1 <- matrix(0,nrow=dim(mu.1)[1],ncol=dim(mu.1)[2])
PP.1 <- rep(0,length(P.1))
for (t in 1:500) {
	print (t)
	mu.1 <- getMu(X.1, factor, Z.1, P.1)
	mu.0 <- getMu(X.0, 1, Z.0, P.0)
	Z.1 <- getZ(X.1, P.1, mu.1, ztable)
	Z.0 <- getZ(X.0, P.0, mu.0, ztable)
	P.1 <- getP(Z.1)
	P.0 <- getP(Z.0)
	mumu.1 <- mumu.1 + mu.1
	PP.1 <- PP.1 + P.1
}
mumu.1 <- mumu.1/500
PP.1 <- PP.1/500

#test stability 2
mumu.2 <- matrix(0,nrow=dim(mu.1)[1],ncol=dim(mu.1)[2])
PP.2 <- rep(0,length(P.1))
prob <- 0
for (t in 1:500) {
	print (t)
	mu.1 <- getMu(X.1, factor, Z.1, P.1)
	mu.0 <- getMu(X.0, 1, Z.0, P.0)
	Z.1 <- getZ(X.1, P.1, mu.1, ztable)
	Z.0 <- getZ(X.0, P.0, mu.0, ztable)
	P.1 <- getP(Z.1)
	P.0 <- getP(Z.0)
	mumu.2 <- mumu.2 + mu.1
	PP.2 <- PP.2 + P.1
	prob <- prob + getProb(test, mu.1, mu.0, Z.1, Z.0)
}
mumu.2 <- mumu.2/500
PP.2 <- PP.2/500
prob <- prob/500