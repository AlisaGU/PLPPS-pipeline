sstselfEval<-function(type,filename){
  # pdf(paste(filename,".pdf",sep=""),width=6,height=6)
  # sink(paste(filename,".txt",sep=""))
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
                     stage=as.numeric(phe[,3]),
                     grade=as.numeric(phe[,4]),
                     score=factor((score[-po4])>(medianScore)))
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
                     score=factor((score.test)>(medianScore)))
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
                     score=factor((score.test.independ)>(medianScore)))
  }else if(type=="global"){
    phe1<-c()
    name1<-c()
    for(i in 1:15){
      x<-clinical[[i]]
      phe1<-rbind(phe1,cbind(x$days_to_death,
                           x$vital_status,
                           x$tumorstage,x$grade))
      name1<-c(name1,rownames(x))
    }
    phe1[,2]<-ifelse(phe1[,2]=="0"|phe1[,2]=="living",0,1)
    
    i=16
    x2<-clinical[[i]]
    phe2<-cbind(x2$donor_survival_time,x2$donor_vital_status,
                x2$tumour_stage,x2$tumour_grade)
    phe2[,2]<-ifelse(phe2[,2]=="1",0,1)
    phe<-rbind(phe1,phe2)
    rownames(phe)<-c(name1,rownames(clinical[[16]]))
    colnames(phe)<-c("time","status","stage","grade")
    data<-data.frame(time=as.numeric(phe[,1]),
                     status=as.numeric(phe[,2]),
                     stage=as.numeric(phe[,3]),
                     grade=as.numeric(phe[,4]),
                     score=factor((-score.all)>(-medianScore)))
  }else if(type=="TestALL"){
    phe.test1<-c()
    name.test1<-c()
    for(i in 1:length(testPo)){
      x<-clinical[[testPo[i]]]
      phe.test1<-rbind(phe.test1,cbind(x$days_to_death,
                           x$vital_status,
                           x$tumorstage,x$grade))
      name.test1<-c(name.test1,rownames(x))
    }
    phe.test1[,2]<-ifelse(phe.test1[,2]=="living",0,1)
    colnames(phe.test1)<-c("time","status","stage","grade")
    rownames(phe.test1)<-name.test1
    
    
    x1<-clinical[[independPo[1]]]
    phe1.test2<-cbind(x1$days_to_death,x1$vital_status,
                x1$tumorstage,x1$grade)
    x2<-clinical[[independPo[[2]]]]
    phe2.test2<-cbind(x2$donor_survival_time,x2$donor_vital_status,
                x2$tumour_stage,x2$tumour_grade)
    phe2.test2[,2]<-ifelse(phe2.test2[,2]==2,1,0)
    phe.test2<-rbind(phe1.test2,phe2.test2)
    rownames(phe.test2)<-c(rownames(x1),rownames(x2))
    colnames(phe.test2)<-c("time","status","stage","grade")
    data<-data.frame(time=c(as.numeric(phe.test1[,1]),as.numeric(phe.test2[,1])),
                     status=c(as.numeric(phe.test1[,2]),as.numeric(phe.test2[,2])),
                     stage=c(as.numeric(phe.test1[,3]),as.numeric(phe.test2[,3])),
                     grade=c(as.numeric(phe.test1[,4]),as.numeric(phe.test2[,4])),
                     score=factor((-c(score.test,score.test.independ))>(-medianScore)))
  }
  #all
  all<-c(as.numeric(table(data$score)),
         summary(coxph(Surv(as.numeric(time),status)~ score,data))$coefficients[-4])
  #T3
  po<-which(data$stage==3)
  data1<-data[po,]
  T3<-c(as.numeric(table(data1$score)),
        summary(coxph(Surv(as.numeric(time),status)~ score,data1))$coefficients[-4])
  data1$score<-ifelse(data1$score==T,"High Risk","Low Risk")
  fit<-survfit(Surv(time,status)~score,data=data1)
  diff<-survdiff(Surv(time,status)~score,data=data1)
  print(fit)
  print(diff)
  plotSurvivalCurve(fit,data1,"T3")
  #T4
  po<-which(data$stage==4)
  data1<-data[po,]
  T4<-c(as.numeric(table(data1$score)),
        summary(coxph(Surv(as.numeric(time),status)~ score,data1))$coefficients[-4])
  data1$score<-ifelse(data1$score==T,"High Risk","Low Risk")
  fit<-survfit(Surv(time,status)~score,data=data1)
  diff<-survdiff(Surv(time,status)~score,data=data1)
  print(fit)
  print(diff)
  plotSurvivalCurve(fit,data1,"T4")
  #G2
  po<-which(data$grade==2)
  data1<-data[po,]
  G2<-c(as.numeric(table(data1$score)),
        summary(coxph(Surv(as.numeric(time),status)~ score,data1))$coefficients[-4])
  data1$score<-ifelse(data1$score==T,"High Risk","Low Risk")
  fit<-survfit(Surv(time,status)~score,data=data1)
  diff<-survdiff(Surv(time,status)~score,data=data1)
  print(fit)
  print(diff)
  plotSurvivalCurve(fit,data1,"G2")
  #G3
  po<-which(data$grade==3)
  data1<-data[po,]
  G3<-c(as.numeric(table(data1$score)),
        summary(coxph(Surv(as.numeric(time),status)~ score,data1))$coefficients[-4])
  data1$score<-ifelse(data1$score==T,"High Risk","Low Risk")
  fit<-survfit(Surv(time,status)~score,data=data1)
  diff<-survdiff(Surv(time,status)~score,data=data1)
  print(fit)
  print(diff)
  plotSurvivalCurve(fit,data1,"G3")
  
  result<-rbind(all,T3,T4,G2,G3)
  colnames(result)<-c("Low Risk","High Risk","lnHR","HR","selnHR","p value")
  # sink()
  # dev.off()
  return(result)
}
