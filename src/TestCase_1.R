#Test case 1
theta_new <- getTheta_MC(alpha = alpha, nVal = nVal)
theta_1 <- theta_new$theta_1
appbest_pept <- matrix(NA, nrow = 1, ncol = 38)
 for(i in 3:1) {
    col_theta <- theta_1[,paste('L',i,sep="")]
    s_class <- which(col_theta == max(col_theta))[1]
    appbest_pept[1,20-i] <- s_class
}
for(i in 1:9) {
    col_theta <- theta_1[,paste('R',i,sep="")]
    s_class <- which(col_theta == max(col_theta))
    appbest_pept[1,19+i] <- s_class
}
prior.positive <- c(10**-(4:1),1/2)
appbest_pred <- c()
for (i in 1:length(prior.positive)) {
    appbest_predict <- c()
    for(j in 1:1000) {
        theta <- getTheta_MC(alpha = alpha, nVal = nVal)
        appbest_predict <- c(appbest_predict,
            NB_predict(appbest_pept, theta, prior.positive = prior.positive[i]))
    }
	appbest_pred <- c(appbest_pred, mean(appbest_predict))
}
plot(log(prior.positive),appbest_pred)
