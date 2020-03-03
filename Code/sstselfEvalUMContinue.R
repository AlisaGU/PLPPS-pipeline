sstselfEvalUMContinue<-function(type){
  meta.training.unicox.score<-c()
  meta.training.unicox.stage<-c()
  meta.training.unicox.grade<-c()
  meta.training.multicox<-c()
  data<-c()
  if(type=="trainPo"){
    phe<-c()
    name<-c()
    for(i in 1:length(trainPo)){
      x<-clinical[[trainPo[i]]]
      phe<-rbind(phe,cbind(x$days_to_death,
                           x$vital_status,
                           x$tumorstage,x$grade))
      name<-c(name,rownames(x))
    }
    phe[,2]<-ifelse(phe[,2]=="living",0,1)
    colnames(phe)<-c("time","status","stage","grade")
    rownames(phe)<-name
    po4<-which(phe[,4]==4)
    phe<-phe[-po4,]
    data<-data.frame(time=as.numeric(phe[,1]),
                     status=as.numeric(phe[,2]),
                     stage=factor(as.numeric(phe[,3])),
                     grade=factor(as.numeric(phe[,4])),
                     score=(scale(-score[-po4])))
  }else if(type=="testPo"){
    phe<-c()
    name<-c()
    for(i in 1:length(testPo)){
      x<-clinical[[testPo[i]]]
      phe<-rbind(phe,cbind(x$days_to_death,
                           x$vital_status,
                           x$tumorstage,x$grade))
      name<-c(name,rownames(x))
    }
    phe[,2]<-ifelse(phe[,2]=="living",0,1)
    colnames(phe)<-c("time","status","stage","grade")
    rownames(phe)<-name
    data<-data.frame(time=as.numeric(phe[,1]),
                     status=as.numeric(phe[,2]),
                     stage=as.numeric(phe[,3]),
                     grade=as.numeric(phe[,4]),
                     score=scale((-score.test)))
  }else if(type=="independent"){
    x1<-clinical[[independPo[1]]]
    phe1<-cbind(x1$days_to_death,x1$vital_status,
                x1$tumorstage,x1$grade)
    x2<-clinical[[independPo[[2]]]]
    phe2<-cbind(x2$donor_survival_time,x2$donor_vital_status,
                x2$tumour_stage,x2$tumour_grade)
    phe2[,2]<-ifelse(phe2[,2]==2,1,0)
    phe<-rbind(phe1,phe2)
    rownames(phe)<-c(rownames(x1),rownames(x2))
    colnames(phe)<-c("time","status","stage","grade")
    data<-data.frame(time=as.numeric(phe[,1]),
                     status=as.numeric(phe[,2]),
                     stage=as.numeric(phe[,3]),
                     grade=as.numeric(phe[,4]),
                     score=scale((-score.test.independ)))
  }
  meta.training.unicox.score<-summary(coxph(Surv(as.numeric(time),status)
                                            ~ score, data))$coefficients
  meta.training.unicox.stage<-summary(coxph(Surv(as.numeric(time),status)
                                            ~ stage, data))$coefficients
  meta.training.unicox.grade<-summary(coxph(Surv(as.numeric(time),status)
                                            ~ grade, data))$coefficients
  meta.training.multicox<-summary(coxph(Surv(as.numeric(time),status)
                                        ~score+stage+grade,data))$coefficients
  unicox<-rbind(meta.training.unicox.score,meta.training.unicox.stage,
                meta.training.unicox.grade)
  rownames(unicox)<-c("HighRisk vs LowRisk","T4 vs T3","G3 vs G2")
  rownames(meta.training.multicox)<-c("HighRisk vs LowRisk","T4 vs T3","G3 vs G2")
  result<-list(unicox=unicox,
               multicox=meta.training.multicox)
  return(result)
}
