data<-read.csv('../../data/binaryData_v2.csv',header=TRUE)
attach(data)
len = c();
for (i in 1:length(AcpH)) {
  len[i]=nchar(as.character(nterm[i]))+nchar(as.character(cterm[i]))+1
}
hitlength0 = len[AcpH==1]
pdf('HitLength0.pdf',width=5,height=3)
hist(hitlength0,xlab=NULL,ylab='# Hits',main=NULL, breaks=1:40, ylim=c(0,45), col="black", cex.axis=1.3, cex.lab=1.3)
dev.off()
detach(data)

data<-read.csv('../../data/newData.csv',header=TRUE)
attach(data)
len = c();
for (i in 1:length(AcpH)) {
  len[i]=nchar(as.character(nterm[i]))+nchar(as.character(cterm[i]))+1
}
hitlength1=c(len[AcpH==1],hitlength0) # newData doesn't include binaryData_v2
pdf('HitLength1.pdf',width=5,height=3)
hist(hitlength1,xlab=NULL,ylab='# Hits',main=NULL, breaks=1:40, ylim=c(0,45), col="black", cex.axis=1.3, cex.lab=1.3)
dev.off()
detach(data)

data<-read.csv('../../data/newData#2.csv',header=TRUE)
attach(data)
len = c();
for (i in 1:length(AcpH)) {
  len[i]=nchar(as.character(nterm[i]))+nchar(as.character(cterm[i]))+1
}
hitlength2=c(len[AcpH==1],hitlength0) # newData2 doesn't include binaryData_v2
pdf('HitLength2.pdf',width=5,height=3)
hist(hitlength2,xlab='Peptide Length',ylab='# Hits',main=NULL, breaks=1:40, ylim=c(0,45), col="black", cex.axis=1.3, cex.lab=1.3)
dev.off()
detach(data)
