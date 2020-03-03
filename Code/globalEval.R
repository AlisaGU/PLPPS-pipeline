globalEval<-function(score,setName,phe,timeCol,eventCol){
  time<- as.numeric(phe[,timeCol])
  event<- as.character(phe[,eventCol])
  event<-ifelse(phe[,eventCol]=="Dead",1,0)
  
  phe<-data.frame(time=time,event=event)
  #score<-(-score)
  alone_score<-summary(coxph(Surv(phe[,timeCol],phe[,eventCol] == 1)
                             ~ score, phe))$coefficients
  roc.ci<-ROC_CI_global(score,phe,setName,timeCol,eventCol)
  HR<-alone_score[2]
  lnHR<-alone_score[1]
  selnHR<-alone_score[3]
  low.HR<-exp(lnHR-1.96*selnHR)
  high.HR<-exp(lnHR+1.96*selnHR)
  p.HR<-alone_score[5]
  oneYearSurvival<-roc.ci[1]
  twoYearSurvival<-roc.ci[2]
  threeYearSurvival<-roc.ci[3]
  fourYearSurvival<-roc.ci[4]
  fiveYearSurvival<-roc.ci[5]
  C_index<-roc.ci[6]
  se_C_index<-roc.ci[7]
  low_C_index<-C_index-1.96*se_C_index
  high_C_index<-C_index+1.96*se_C_index
  result<-structure(c(HR,lnHR,selnHR,low.HR,high.HR,p.HR,oneYearSurvival,twoYearSurvival,
                      threeYearSurvival,fourYearSurvival,fiveYearSurvival,
                      C_index,se_C_index,low_C_index,high_C_index),
                    names=c("HR","lnHR","selnHR","low.HR",
                            "high.HR","p.HR","oneYearSurvival","twoYearSurvival",
                            "threeYearSurvival","fourYearSurvival","fiveYearSurvival","C_index",
                            "se_C_index","low_C_index","high_C_index"))
  return(result)
}