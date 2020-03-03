ROC3_5_CI<-function(score,phe,timeCol,eventCol){
  time<- as.numeric(phe[,timeCol])
  event<- as.character(phe[,eventCol])
  event<-ifelse(phe[,eventCol]=="Dead",1,0)
  
  phe<-data.frame(time=time,event=event)
  a<-timeROC(T=phe$time,
             delta=phe$event,marker=score,
             cause=1,weighting = "marginal",
             times=c(365.25*3,365.25*5),iid = T)
  b<-c(a$AUC[1],a$inference$vect_sd_1[1],a$AUC[2],a$inference$vect_sd_1[2])
  names(b)<-c("aucOf3","se.aucof3","aucOf5","se.aucof5")
  return(b)
}