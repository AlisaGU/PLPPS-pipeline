selfEval<-function(score,phe,setName,sStatus,gStatus,aStatus,sCol,gCol,aCol,timeCol,eventCol){
  phe[,timeCol]<- as.numeric(phe[,timeCol])
  phe[,eventCol]<- as.character(phe[,eventCol])
  if(any("deceased"%in%phe[,eventCol])){
    phe[,eventCol]<-ifelse(phe[,eventCol]=="deceased",1,0)
  }
  #score<-(-score)
  #uniCox
  alone_score<-summary(coxph(Surv(phe[,timeCol],phe[,eventCol] == 1)
                                        ~ score, phe))$coefficients
  alone_tumorstage<-c()
  alone_age<-c()
  alone_grade<-c()
  if("III"%in%phe[,sCol]){
    phe[,sCol]<-ifelse(phe[,sCol]=="III",3,4)
  }
  if(sStatus){
    alone_tumorstage<-summary(coxph(Surv(phe[,timeCol],phe[,eventCol] == 1)~
                                      phe[,sCol],phe))$coefficients
  }
  if(gStatus){
    alone_grade<-summary(coxph(Surv(phe[,timeCol],phe[,eventCol] == 1)~
                                      phe[,gCol],phe))$coefficients
  }
  if(aStatus){
    alone_age<-summary(coxph(Surv(phe[,timeCol],phe[,eventCol] == 1)~
                                      phe[,aCol],phe))$coefficients
  }
  
  ##ROC-CI
  roc.ci<-ROC_CI(score,phe,setName,timeCol,eventCol)
  
  ##multiCox
  multiFactors<-c()
  if(sStatus){
    if(gStatus&aStatus){#8,11,12,14,15,16
      multiFactors <- summary(coxph(Surv(phe[,timeCol],phe[,eventCol] == 1)
                                    ~ score+phe[,sCol]+phe[,gCol]+
                                      phe[,aCol], phe))$coefficients
    }
    if(!gStatus&aStatus){#1,7
      multiFactors <- summary(coxph(Surv(phe[,timeCol],phe[,eventCol] == 1)
                                    ~ score+phe[,sCol]+phe[,aCol], phe))$coefficients
    }
    if(gStatus& !aStatus){#3,4,6,9,10,13
      multiFactors <- summary(coxph(Surv(phe[,timeCol],phe[,eventCol] == 1)
                                    ~ score+phe[,sCol]+phe[,gCol], phe))$coefficients
    }
  }else{#2,5
    if(gStatus){#2
      multiFactors <- summary(coxph(Surv(phe[,timeCol],phe[,eventCol] == 1)
                                    ~ score+phe[,gCol]+phe[,aCol], phe))$coefficients
    }else{#5
      multiFactors <- summary(coxph(Surv(phe[,timeCol],phe[,eventCol] == 1)
                                    ~ score, phe))$coefficients
    }
  }
  return(list(multiCox=multiFactors,
              uniCox=rbind(alone_score,alone_tumorstage,alone_age,alone_grade),
              roc_ci=roc.ci))
}