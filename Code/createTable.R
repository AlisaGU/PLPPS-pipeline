createTable<-function(data){
  #·Ö±ðÎªlnHR,HR,selnHR,z,p
  HR<-data$`exp(coef)`
  lnHR<-data$coef
  selnHR<-data$`se(coef)`
  low.HR<-exp(lnHR-1.96*selnHR)
  high.HR<-exp(lnHR+1.96*selnHR)
  a<-data.frame(Variable=rownames(data),
                format(HR,digits=2),
                paste(format(low.HR,digits=2),
                      format(high.HR,digits=3),sep="-"),
                p=format(data$`Pr(>|z|)`,digits=2,scientific = 2))
  colnames(a)<-c("Variable","HR","95% CI of HR","P value")
  return(a)
}