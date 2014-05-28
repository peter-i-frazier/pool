data <- read.csv('../../src/compare_maxp_vs_naive/forplot.csv')
# The first row of data is maxp, the second row is ranking by probability of a
# hit, the third row is mutate.

pdf('maxp_horserace.pdf',width=4,height=3)
plot(x=1:100,y=data[2,2:101],xlab='# Peptides Tested',ylab='Probability of a short hit',type='l',ylim=c(0,1),cex=1.5,cex.lab=1.5,lwd=2,col="green",cex.axis=1.5,axes=FALSE) # ranking
axis(1,at=c(0,20,40,60,80,100),cex.axis=1.5)
axis(2,at=c(0.0,0.5,1.0),cex.axis=1.5)
lines(x=1:100,y=data[1,2:101],col="red",lwd=2) # maxp
lines(x=1:100,y=data[3,2:101],col="blue",lwd=2) # mutate
# leg.txt <- c('Basic Ranking', 'POOL') 
dev.off()


# Display the same plot but without maxp
pdf('horserace_without_maxp.pdf',width=4,height=3)
plot(x=1:100,y=data[2,2:101],xlab='# Peptides Tested',ylab='Probability of a short hit',type='l',ylim=c(0,1),cex=1.5,cex.lab=1.5,lwd=2,col="green",cex.axis=1.5,axes=FALSE) # ranking
axis(1,at=c(0,20,40,60,80,100),cex.axis=1.5)
axis(2,at=c(0.0,0.5,1.0),cex.axis=1.5)
#lines(x=1:100,y=data[1,2:101],col="red",lwd=2) # maxp
lines(x=1:100,y=data[3,2:101],col="blue",lwd=2) # mutate
# leg.txt <- c('Basic Ranking', 'POOL') 
dev.off()

