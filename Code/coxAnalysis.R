coxAnalysis<-function(Mat,phe,stageStatus,gradeStatus,ageStatus,expNAstatus){
  a<-c()
  d<-c()
  if(expNAstatus==1){#15,14
    for (j in 1:nrow(Mat)) {
      expj <- Mat[j,];
      if(length(expj[!is.na(expj)])!=0){
        gene.j <- coxph(Surv(days_to_death,vital_status == "deceased")
                        ~ expj[!is.na(expj)]+tumorstage+grade+age_at_initial_pathologic_diagnosis, phe[!is.na(expj),])
        b<-summary(gene.j)
        a<-rbind(a,c(b$coefficients[1,5],b$conf.int[1,c(1,3,4)]))
        d<-rbind(d,b$coefficients[1,1:3])
      }
    }
  }else if(expNAstatus==0&stageStatus==T){#1:15[-C(2,5)]
    if(gradeStatus & ageStatus) {
      for (j in 1:nrow(Mat)) {
        gene.j <- coxph(Surv(days_to_death,vital_status == "deceased")
                        ~ Mat[j,]+tumorstage+grade+age_at_initial_pathologic_diagnosis, phe)
        b<-summary(gene.j)
        a<-rbind(a,c(b$coefficients[1,5],b$conf.int[1,c(1,3,4)]))
        d<-rbind(d,b$coefficients[1,1:3])
      }
    }
    if(!gradeStatus & ageStatus) {
      for (j in 1:nrow(Mat)) {
        gene.j <- coxph(Surv(days_to_death,vital_status == "deceased")
                        ~ Mat[j,]+tumorstage+age_at_initial_pathologic_diagnosis, phe)
        b<-summary(gene.j)
        a<-rbind(a,c(b$coefficients[1,5],b$conf.int[1,c(1,3,4)]))
        d<-rbind(d,b$coefficients[1,1:3])
      }
      
    }
    if (gradeStatus & !ageStatus) {
      for (j in 1:nrow(Mat)) {
        gene.j <- coxph(Surv(days_to_death,vital_status == "deceased")
                        ~ Mat[j,]+tumorstage+grade, phe)
        b<-summary(gene.j)
        a<-rbind(a,c(b$coefficients[1,5],b$conf.int[1,c(1,3,4)]))
        d<-rbind(d,b$coefficients[1,1:3])
      }
    }
    if (!gradeStatus & !ageStatus) {
      for (j in 1:nrow(Mat)) {
        gene.j <- coxph(Surv(days_to_death,vital_status == "deceased")
                        ~ Mat[j,]+tumorstage, phe)
        b<-summary(gene.j)
        a<-rbind(a,c(b$coefficients[1,5],b$conf.int[1,c(1,3,4)]))
        d<-rbind(d,b$coefficients[1,1:3])
      }
    }
  }else if(expNAstatus==0&stageStatus==F){#2,5
    if(gradeStatus==T){#2
      for (j in 1:nrow(Mat)) {
        gene.j <- coxph(Surv(days_to_death,vital_status == "deceased")
                        ~ Mat[j,]+grade+age_at_initial_pathologic_diagnosis, phe)
        b<-summary(gene.j)
        a<-rbind(a,c(b$coefficients[1,5],b$conf.int[1,c(1,3,4)]))
        d<-rbind(d,b$coefficients[1,1:3])
      }
    }else{#5
      for (j in 1:nrow(Mat)) {
        gene.j <- coxph(Surv(days_to_death,vital_status == "deceased")
                        ~ Mat[j,], phe)
        b<-summary(gene.j)
        a<-rbind(a,c(b$coefficients[1,5],b$conf.int[1,c(1,3,4)]))
        d<-rbind(d,b$coefficients[1,1:3])
      }
    }
  }
  po<-apply(Mat,1,function(x){!all(is.na(x))})
  rownames(a)<-rownames(Mat)[po]
  rownames(d)<-rownames(Mat)[po]
  colnames(a)<-c("cox.p","HR","95%CI.low","95%CI.high")
  colnames(d)<-c("lnHR","HR","se(lnHR)")
  result<-list(a,d)
  return(result)
}
