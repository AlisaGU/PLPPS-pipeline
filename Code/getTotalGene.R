getTotalGene<-function(x){
  #x是一个列表
  re<-rownames(x[[1]])
  for(i in 2:length(x)){
    re<-union(re,rownames(x[[i]]))
  }
  return(re)
}