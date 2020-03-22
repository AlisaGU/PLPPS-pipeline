MulticoxForestData<-function(Data){
  HR<-sapply(Data,function(x){x[[1]][1,2]})
  lnHR<-sapply(Data,function(x){x[[1]][1,1]})
  selnHR<-sapply(Data,function(x){x[[1]][1,3]})
  low.HR<-exp(lnHR-1.96*selnHR)
  high.HR<-exp(lnHR+1.96*selnHR)
  C_index<-sapply(Data,function(x){x[[3]][3]})
  se_C_index<-sapply(Data,function(x){x[[3]][4]})
  low_c_index<-C_index-1.96*se_C_index
  
  high_c_index<-C_index+1.96*se_C_index
  HR_C<-data.frame(HR,lnHR,selnHR,low.HR,high.HR,
                   C_index,se_C_index,low_c_index,high_c_index)
  lnHR<-metagen(HR_C[,2],HR_C[,3],sm="lnHR")
  if(lnHR$pval.Q>0.05){
    summary_HR<-c(lnHR$TE.fixed,lnHR$seTE.fixed,lnHR$pval.fixed)
  }else{
    summary_HR<-c(lnHR$TE.random,lnHR$seTE.random,lnHR$pval.random)
  }
  C<-metagen(HR_C[,6],HR_C[,7],sm="C-index",null.effect=0.5)
  if(C$pval.Q>0.05){
    summary_C<-c(C$TE.fixed,C$seTE.fixed,C$pval.fixed)
  }else{
    summary_C<-c(C$TE.random,C$seTE.random,C$pval.random)
  }
  # lnHR<-meta.summaries(HR_C[,2],HR_C[,3],method="fixed",logscale=T)
  # if(lnHR$het[3]>0.05){
  #   summary_HR<-c(lnHR$summary,lnHR$se.summary)
  # }else{
  #   lnHR<-meta.summaries(HR_C[,2],HR_C[,3],method="random",logscale=T)
  #   summary_HR<-c(lnHR$summary,lnHR$se.summary)
  # }
  # #求C-statics的meta结果
  # C<-meta.summaries(HR_C[,6],HR_C[,7],method="fixed",logscale=F)
  # if(C$het[3]>0.05){
  #   summary_C<-c(C$summary,C$se.summary)
  # }else{
  #   C<-meta.summaries(HR_C[,6],HR_C[,7],method="random",logscale=F)
  #   summary_C<-c(C$summary,C$se.summary)
  # }
  ssgsea_colla<-c(summary_HR,summary_C)
  return(list(ssgsea_colla,HR_C))
}