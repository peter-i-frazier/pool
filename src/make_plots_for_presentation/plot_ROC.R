# this is the ROC curve for training set 1
data<-read.csv('../../result/ROC_curve/ROC_DS1_1000_0.25.csv',header=TRUE)
# data contains the true positive rate (TPR) and false positive rate (FPR).

pdf('2013_12_AFOSR_ROC_DS1_small.pdf',width=3,height=3)
plot(data$FPR,data$TPR,type="l",xlab=FALSE,ylab=FALSE,cex=2,xlim=c(0,1),ylim=c(0,1),ann=FALSE,lwd=4,axes=FALSE,frame.plot=TRUE)
axis(1,at=c(0,1),cex.axis=2)
axis(2,at=c(0,1),cex.axis=2)
dev.off()

pdf('2013_12_AFOSR_ROC_DS1_big.pdf',width=5,height=5)
plot(data$FPR,data$TPR,type="l",xlab='False Positive Rate',ylab='True Positive Rate',cex.lab=1.5,xlim=c(0,1),ylim=c(0,1),lwd=2,axes=FALSE,frame.plot=TRUE)
axis(1,at=c(0,1),cex.axis=1.5)
axis(2,at=c(0,1),cex.axis=1.5)
dev.off()


