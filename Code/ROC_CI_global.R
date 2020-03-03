ROC_CI_global<-function(score,phe,setName,timeCol,eventCol){
  alone_roc1<-survivalROC(phe[,timeCol],phe[,eventCol] == 1,
                         (score),predict.time = 365.25*1,method = "KM")
  plotSurvivalRoc(alone_roc1,paste(setName,":method= KM \n Year=1"))
  
  alone_roc2<-survivalROC(phe[,timeCol],phe[,eventCol] == 1,
                         (score),predict.time = 365.25*2,method = "KM")
  plotSurvivalRoc(alone_roc2,paste(setName,":method= KM \n Year=2"))
  
  alone_roc3<-survivalROC(phe[,timeCol],phe[,eventCol] == 1,
                          (score),predict.time = 365.25*3,method = "KM")
  plotSurvivalRoc(alone_roc3,paste(setName,":method= KM \n Year=3"))
  
  alone_roc4<-survivalROC(phe[,timeCol],phe[,eventCol] == 1,
                         (score),predict.time = 365.25*4,method = "KM")
  plotSurvivalRoc(alone_roc4,paste(setName,":method= KM \n Year=4"))
  
  alone_roc5<-survivalROC(phe[,timeCol],phe[,eventCol] == 1,
                          (score),predict.time = 365.25*5,method = "KM")
  plotSurvivalRoc(alone_roc5,paste(setName,":method= KM \n Year=5"))
  # 
  # alone_roc10<-survivalROC(phe[,timeCol],phe[,eventCol] == 1,
  #                         (score),predict.time = 365.25*10,method = "KM")
  # plotSurvivalRoc(alone_roc10,paste(setName,":method= KM \n Year=10"))
  alone_roc_ci<-c(alone_roc1$AUC,alone_roc2$AUC,alone_roc3$AUC,alone_roc4$AUC,alone_roc5$AUC)
  alone_c<-survConcordance(Surv(phe[,timeCol],phe[,eventCol] == 1) ~(score))
  return(c(alone_roc_ci,alone_c$concordance,alone_c$std.err))
}