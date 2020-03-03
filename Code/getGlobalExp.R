getGlobalExp<-function(genes,exp){
  #exp是一个列表
  re<-c()
  for(i in 1:length(exp)){
    genes.miss<-genes[is.na(match(genes,rownames(exp[[i]])))]
    missExp<-matrix(NA,nrow=length(genes.miss),ncol=ncol(exp[[i]]))
    rownames(missExp)<-genes.miss
    colnames(missExp)<-colnames(exp[[i]])
    a<-rbind(exp[[i]],missExp)
    b<-a[match(genes,rownames(a)),]
    re<-cbind(re,b)
  }
  return(re)
}
