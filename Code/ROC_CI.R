ROC_CI<-function(score,phe,setName,timeCol,eventCol){
  alone_roc3<-survivalROC(phe[,timeCol],phe[,eventCol] == 1,
                          (score),predict.time = 365.25*3,method = "KM")
  #plotSurvivalRoc(alone_roc3,paste(setName,":method= KM \n Year=3"))
  alone_roc5<-survivalROC(phe[,timeCol],phe[,eventCol] == 1,
                          (score),predict.time = 365.25*5,method = "KM")
  #plotSurvivalRoc(alone_roc5,paste(setName,":method= KM \n Year=5"))
  alone_roc_ci<-c(alone_roc3$AUC,alone_roc5$AUC)
  alone_c<-survConcordance(Surv(phe[,timeCol],phe[,eventCol] == 1) ~(score))
  return(c(alone_roc_ci,alone_c$concordance,alone_c$std.err))
}