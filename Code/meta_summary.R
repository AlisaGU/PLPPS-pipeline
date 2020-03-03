meta_summary<-function(Metaset){#训练集中满足频率条件的基因的lnHR,HR,se(lnHR)
  gene_meta_summary<-c()
  for(i in 1:length(Metaset)){
    logHR<-Metaset[[i]][,1]
    selogHR<-Metaset[[i]][,3]
    resultMeta<-metagen(logHR,selogHR,sm="HR")
    summary.value<-c(resultMeta$Q, resultMeta$df.Q,1-pchisq(resultMeta$Q, resultMeta$df.Q),
                     resultMeta$I2,exp(resultMeta$TE.fixed),exp(resultMeta$upper.fixed),
                     exp(resultMeta$lower.fixed),2*pnorm(abs(resultMeta$zval.fixed),lower.tail=FALSE),
                     exp(resultMeta$TE.random),exp(resultMeta$upper.random),
                     exp(resultMeta$lower.random),2*pnorm(abs(resultMeta$zval.random),lower.tail=FALSE))
    gene_meta_summary<-rbind(gene_meta_summary,summary.value)
  }
  rownames(gene_meta_summary)<-names(Metaset)
  colnames(gene_meta_summary)<-c("Q","df.Q","p.heterogeneity",
                                 "I2","fixed.HR","fixed.upper",
                                 "fixed.lower","fixed.p","random.HR",
                                 "random.upper","random.lower","random.p")
  return(gene_meta_summary)
}