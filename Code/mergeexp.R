mergeexp<-function(list){
  data<-c()
  for(i in 1:(length(list)-1)){
    data<-merge(list[[i]],list[[i+1]],all=T)
  }
}