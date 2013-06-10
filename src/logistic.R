library(gdata)
rm(list=ls())

AAclass <- read.xls('../data/AAclass.xlsx', header=T, as.is=T)
binary.data <- data.frame(read.xls('../data/binaryData.xlsx', header=T, as.is=T))
ez2.feature <- matrix(0, nrow=dim(binary.data[1]), ncol=6*38)
for (r in 1:dim(binary.data)[1]) {
	sequence <- unlist(strsplit(binary.data[r, 'nterm'],''))
	for (i in 1:length(sequence)) {
		c.no <- (i-1)*6 + AAclass[1,sequence[i]]
		ez2.feature[r, c.no] <- 1
	}
	sequence <- unlist(strsplit(binary.data[r, 'cterm'],''))
	for (i in 1:length(sequence)) {
		c.no <- (i-1+19)*6 + AAclass[1,sequence[i]]
		ez2.feature[r, c.no] <- 1
	}
}
ez2 <- data.frame(cbind(binary.data[,2],ez2.feature))
names(ez2)[1] <- 'y'

MAP <- function(x, train.x, train.y) {
	p <- 0.8
	M <- 10000
	sigma0 <- 10
	sigmab <- x[1]
	mu <- x[-1]
	theta.table <- rnorm(M, mean=0, sd=sigma0)
	for (i in 1:dim(train.x)[2]) theta.table <- rbind(theta.table, rnorm(M,mean=mu[i], sd=sigmab))
	X <- cbind(rep(1,dim(train.x)[1]), train.x)
	prob.table <- X %*% theta.table
	prob.table <- prob.table >=0
	Q <- p * prob.table + (1-p) * (1-prob.table)
	v <- train.y %*% log(Q) + (1-train.y) %*% log(1-Q)
	mv <- -mean(v)
}

x <- c(1,rep(1,229))
to.optimize <- function(param) return (MAP(param, ez2.feature, ez2[,'y']))
opt <- optim(par=x, fn=to.optimize, control=list(maxit=2000))

