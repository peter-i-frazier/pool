# Author: Tom Fei
# Created: 06.06.2015
# Plot curvature final

plot.data <- read.csv("curv_final_set.csv")
plot.data <- unlist(unname(plot.data))
#png(filename = "curv_final_set.png")
#pdf("curv_final_set_curv.pdf")
pdf("curv_final_set_hist.pdf")
itr <- 10000
title <- sprintf("Curvature final set, itr = %d
            mean = %.2f, sd = %.2f, 
            min = %.2f, max = %.2f", itr,
                 mean(plot.data), sd(plot.data),
                 min(plot.data), max(plot.data))
#par(mfrow = c(1, 2))
#plot(plot.data, main = title,
#    # type = "l",
#     lwd = 2, xlab = "# of iterations", ylab = "Curvature")
hist(plot.data)
dev.off()


