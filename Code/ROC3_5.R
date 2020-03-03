ROC3_5<-function(score,phe,timeCol,eventCol){
  time<- as.numeric(phe[,timeCol])
  event<- as.character(phe[,eventCol])
  event<-ifelse(phe[,eventCol]=="Dead",1,0)
  
  phe<-data.frame(time=time,event=event)
  alone_roc3<-survivalROC(phe[,timeCol],phe[,eventCol] == 1,
                          (score),predict.time = 365.25*3,method = "KM")
  alone_roc5<-survivalROC(phe[,timeCol],phe[,eventCol] == 1,
                          (score),predict.time = 365.25*5,method = "KM")
  return(c(alone_roc3$AUC,alone_roc5$AUC))
}