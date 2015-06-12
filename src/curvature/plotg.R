# Author: Tom Fei
# Created: 06.06.2015
# Plot curvature growing

plot.data <- read.csv("curv_grow_set.csv")
plot.data <- unlist(unname(plot.data))
png(filename = "curv_grow_set.png")
itr <- 10000
title <- sprintf("Curvature growing set, itr = %d
                 mean = %.2f, sd = %.2f, min = %.2f, max = %.2f", itr,
                 mean(plot.data), sd(plot.data),
                 min(plot.data), max(plot.data))
plot(plot.data, main = title,
     type = "l", lwd = 2, xlab = "# of iterations", ylab = "Curvature")
dev.off()
