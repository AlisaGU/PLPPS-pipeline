chemo.label<-function(x,criterion.left,criterion.right){
  po<-which(as.numeric(x[1])>=criterion.left&as.numeric(x[1])<criterion.right)
  #criterion<-criterion.left[po]
  return(po)
  # if(criterion==0){
  #   return("0")
  # }else{
  #   return(substr(criterion,1,nchar(criterion)-2)) 
  # }
}