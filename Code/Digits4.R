Digits4<-function(x,type="Matrix"){
  if(type=="vector"){
    result<-vector(mode="character",length=length(x))
    for(i in 1:length(x)){
      result[i]<-format(x[i],digits = 4)
    }
    names(result)<-names(x)
    return(result)
  }else{
    result<-matrix(NA,nrow=nrow(x),ncol = ncol(x))
    for(i in 1:nrow(x)){
      for(j in 1:ncol(x)){
        result[i,j]<-format(x[i,j],digits = 4)
      }
    }
    colnames(result)<colnames(x)
    rownames(result)<-rownames(x)
    return(result)
  }
}