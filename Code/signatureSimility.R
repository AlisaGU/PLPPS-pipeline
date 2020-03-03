signatureSimility<-function(sa,sb,globalExp){
  rho<-matrix(NA,nrow=length(sa),ncol=length(sb))
  for(i in 1:length(sa)){
    for(j in 1:length(sb)){
      ea.i<-as.numeric(globalExp[which(sa[i]==rownames(globalExp)),])
      na.i<-which(!is.na(ea.i))
      eb.j<-as.numeric(globalExp[which(sb[j]==rownames(globalExp)),])
      nb.j<-which(!is.na(eb.j))
      ea.i<-ea.i[intersect(na.i,nb.j)]
      eb.j<-eb.j[intersect(na.i,nb.j)]
      rho[i,j]<-cor(ea.i,eb.j)
    }
  }
  acor = as.numeric(abs(rho))
  cutoff = (exp(2*1.96/sqrt(nrow(globalExp)-3))-1)/(exp(2*1.96/sqrt(nrow(globalExp)-3))+1)
  Cor = sum(acor[!is.na(acor)]>cutoff)/prod(c(length(sa), length(sb)))
  return(Cor)
}