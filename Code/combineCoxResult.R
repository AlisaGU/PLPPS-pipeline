combineCoxResult<-function(x,Xname){#eg: selfEvalResult[[1]] 
  multi.lnHR<-x[[1]][,1]
  multi.HR<-x[[1]][,2]
  multi.selnHR<-x[[1]][,3]
  multi.p<-x[[1]][,5]
  multi.low<-exp(multi.lnHR-1.96*multi.selnHR)
  multi.high<-exp(multi.lnHR+1.96*multi.selnHR)
  
  uni.lnHR<-x[[2]][,1]
  uni.HR<-x[[2]][,2]
  uni.selnHR<-x[[2]][,3]
  uni.p<-x[[2]][,5]
  uni.low<-exp(uni.lnHR-1.96*uni.selnHR)
  uni.high<-exp(uni.lnHR+1.96*uni.selnHR)
  
  re<-data.frame(format(multi.HR,digits=4),
                 paste(format(multi.low,digits=4),
                       format(multi.high,digits=4),sep="-"),
                 format(multi.p,scientific =T,digits = 3),
                 format(uni.HR,digits=4),
                 paste(format(uni.low,digits=4),
                       format(uni.high,digits=4),sep="-"),
                 format(uni.p,scientific =T,digits = 3))
  
  rownames(re)[grepl("score",rownames(re))]<-"Score"
  rownames(re)[grepl("sCol",rownames(re))]<-"Stage"
  rownames(re)[grepl("gCol",rownames(re))]<-"Grade"
  rownames(re)[grepl("aCol",rownames(re))]<-"Age"
  if(nrow(re)==1){
    rownames(re)<-paste(Xname,"Score",sep="%")
  }else{
    rownames(re)[1]<-paste(Xname,rownames(re)[1],sep="%")
  }
  re1<-data.frame(rownames(re),re)
  names(re1)<-c("Variables","Multi:HR","Multi:95%CI of HR","Multi:P value",
               "Uni:HR","Uni:95%CI of HR","Uni:P value")
  return(re1)
}