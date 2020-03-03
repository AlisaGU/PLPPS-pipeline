getHeatmapExp<-function(gene,expList){
  Exp<-c()
  for(i in 1:length(expList)){
    Exp<-cbind(Exp,expList[[i]][match(gene,rownames(expList[[i]])),])
  }
  rownames(Exp)<-gene
  colnames(Exp)<-unlist(sapply(expList,colnames))
  return(Exp)
}