ROC_plot <- function(prob, Y, filename, ptitle) 
{
	FPR <- rep(-1, length(prob))
	TPR <- rep(-1, length(prob))
	thresholds <- sort(prob)
	for( i in 1:length(prob) ) {
		threshold <- thresholds[i]
		label <- rep(0, length(prob))
		for( j in 1:length(prob) ){
			if(prob[j] >= threshold) {
				label[j] <- 1 }
		}
		FPR[i] <- sum((Y==0)&(label==1))/sum(Y==0)
		TPR[i] <- sum((Y==1)&(label==1))/sum(Y==1)
	}
	pdf(filename)
	plot(x = FPR, y = TPR, type = 'l', xlim = c(0,1), ylim = c(0,1), xlab = "false positive rate", ylab = "true positive rate")
	title(main = ptitle)
	dev.off()
}