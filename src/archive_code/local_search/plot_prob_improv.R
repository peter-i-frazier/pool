#!/usr/bin/env Rscript

# Plot PI from file
prob.improv.S <- read.csv("prob_improv_S_data.csv", header = F)
png(filename = "plot.png")
plot(prob.improv.S$V1, main = "Local Search", type = "l",
     lwd = 2, xlab = "# of iterations", ylab = "PI(S)")
dev.off()
