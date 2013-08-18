for (t in 1:1000) {
	print (t)
	mu.1 <- getMu(X.1, factor, Z.1, P.1)
	mu.0 <- getMu(X.0, 1, Z.0, P.0)
	Z.1 <- getZ(X.1, P.1, mu.1, ztable)
	Z.0 <- getZ(X.0, P.0, mu.0, ztable)
	P.1 <- getP(Z.1)
	P.0 <- getP(Z.0)
}
prob <- getProb(test, mu.1, mu.0, Z.1, Z.0)