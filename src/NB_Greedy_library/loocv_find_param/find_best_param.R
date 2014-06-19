a<-unlist(read.table('auc',sep='\n'))
print (order(a,decreasing=T)[1])
