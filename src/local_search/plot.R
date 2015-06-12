# Author: Tom Fei
# Created: 06.06.2015
# Plot PI(S) from file

plot.data <- read.csv("prob_improv_S_10000itr.csv")
plot.data <- unlist(unname(plot.data))
png(filename = "plot_pi_10000itr.png")
plot(plot.data, main = "Local Search", type = "l",
     lwd = 2, xlab = "# of iterations", ylab = "PI(S)")
dev.off()
