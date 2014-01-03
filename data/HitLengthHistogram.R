data<-read.csv('../../data/binaryData_v2.csv',header=TRUE)
attach(data)
len = c();
for (i in 1:length(AcpH)) {
  len[i]=nchar(as.character(nterm[i]))+nchar(as.character(cterm[i]))+1
}
pdf('HitLength0.pdf',width=4,height=3)
hist(len[AcpH==1],xlab='Peptide Length',ylab='# Hits',main=NULL, breaks=1:40)
dev.off()
detach(data)



data<-read.csv('../../data/newData.csv',header=TRUE)
attach(data)
len = c();
for (i in 1:length(AcpH)) {
  len[i]=nchar(as.character(nterm[i]))+nchar(as.character(cterm[i]))+1
}
pdf('HitLength1.pdf',width=4,height=3)
hist(len[AcpH==1],xlab='Peptide Length',ylab='# Hits',main=NULL, breaks=1:40)
dev.off()
detach(data)

data<-read.csv('../../data/newData#2.csv',header=TRUE)
attach(data)
len = c();
for (i in 1:length(AcpH)) {
  len[i]=nchar(as.character(nterm[i]))+nchar(as.character(cterm[i]))+1
}
pdf('HitLength2.pdf',width=4,height=3)
hist(len[AcpH==1],xlab='Peptide Length',ylab='# Hits',main=NULL, breaks=1:40)
dev.off()
detach(data)
