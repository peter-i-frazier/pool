# Author: Tom Fei
# Created: 06.06.2015
# Plot PI(S) from file

plot.data <- read.csv("pi_S_500_iter_1000_rerun.csv")
#plot.data <- read.csv("prob_improv_S_1000itr.csv")
plot.data <- unlist(unname(plot.data))
#png(filename = "plot_pi_10000itr.png")
pdf("plot_pi_S_500_iter_1000_rerun.pdf")
plot(plot.data, main = "Local Search, |S| = 500", type = "l",
     lwd = 2, xlab = "# of iterations", ylab = "PI(S)")
dev.off()
