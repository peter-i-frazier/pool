# Author: Tom Fei
# Created: 06.06.2015
# Plot PI(S) from file

plot.data <- read.csv("prob_improv_S_uniq.csv")
plot.data <- unlist(unname(plot.data))
png(filename = "plot_pi_uniq.png")
plot(plot.data, main = "Local Search (Unique)", type = "l",
     lwd = 2, xlab = "# of iterations", ylab = "PI(S)")
dev.off()
