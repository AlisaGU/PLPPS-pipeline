library(survHD)
library(GEOquery)
library(survival)
library(RSQLite)
library(genefilter)
Rcode<-dir("G:\\haodapeng-data\\deleteScatter\\code\\surv")
Rcode<-paste("G:\\haodapeng-data\\deleteScatter\\code\\surv\\",
             Rcode,sep="")
for(i in 1:length(Rcode))
{
  source(Rcode[i]);
  cat(i,"\n")
}

#####t_score#####
t_score<-function(exp,positive,negative)#exp:表达谱矩阵;positive:系数>1的基因;negative:系数<1的基因;
{
  n1<-length(positive);n2<-length(negative)
  apply(exp,2,function(x){
    x1<-x[positive]
    x2<-x[negative]
    -(t.test(x1,x2,alternative = "two.sided")$statistic)
  })
}
#####modelScore#####
modelScore<-function(model_file,input_file)#model_file为模型的名字,input_file为表达谱的名字
{
  model<-model.coefs1[[which(names(model.coefs1)==model_file)]]
  if(model_file=="kernagis12"){
    exp<-datasets_nonscale[[which(names(datasets_nonscale)==input_file)]]
  }else{
    exp<-datasets_scale[[which(names(datasets_scale)==input_file)]]
  }
  if(model_file=="mok09"){
    po<-apply(exp,1,function(x){all(!is.na(x))})
    exp<-exp[po,]
  }
  inters_row<-rownames(exp)[rownames(exp)%in%names(model)]
  inters_col<-colnames(exp)
  
  if(model_file!="denkert09"){
    if(length(inters_row)==1){
      exp<-matrix(exp[inters_row,],
                  nrow=length(inters_row))
    }else{
      exp<-exp[inters_row,]
    }
  }else{
    exp<-2^(matrix(exp[inters_row,],nrow=length(inters_row)))
  }
  rownames(exp)<-inters_row
  colnames(exp)<-inters_col
  if(nrow(exp)==1){
    exp<-matrix(exp[,!is.na(exp)],nrow=1)
    rownames(exp)<-inters_row
    colnames(exp)<-inters_col[-92]
  }
  if(model_file=="ours116"){
    positive<-names(model[model>0])[names(model[model>0])%in%rownames(exp)]
    negative<-names(model[model<0])[names(model[model<0])%in%rownames(exp)]
    model_score<-t_score(exp,positive,negative)
  }else{
    methodModel<-model_method[[which(names(model_method)==model_file)]]
    coef<-structure(methodModel@coefficients[rownames(exp)],
                    .names=rownames(exp))
    methodModel1<-initialize(methodModel,coefficients=coef)
    model_score<-predict(methodModel1, newdata=t(exp), type="lp")@lp
  }
  nonNA_po<-which(!is.na(model_score))
  model_score<-structure(model_score[nonNA_po],names=colnames(exp)[nonNA_po])
  return(model_score)
  
}

#####modelMetaStatics#####
modelMetaStatics<-function(model_file,input_file){
  model_score<-scale(modelScore(model_file,input_file))
  clinical_model<-clinical[[which(input_file==names(clinical))]]
  a<-c()
  if(input_file!="Patch et al."){
    nonNA_po<-match(rownames(model_score),rownames(clinical_model))
    a<-cbind(clinical_model$days_to_death[nonNA_po],
             clinical_model$vital_status[nonNA_po],model_score)
    a[,2]<-ifelse(a[,2]=="living",0,1)
  }else{
    nonNA_po<-match(rownames(model_score),rownames(clinical_model))
    a<-cbind(clinical_model$donor_survival_time[nonNA_po],
             clinical_model$donor_vital_status[nonNA_po],model_score)
    a[,2]<-ifelse(a[,2]=="alive",0,1)
  }
  rownames(a)<-names(model_score)
  colnames(a)<-c("time","status","model_score")
  a<-as.data.frame(apply(a,2,as.numeric))
  bb<-survConcordance(Surv(time,status)~model_score,a)#一致性
  cc<-coxph(Surv(time,status)~model_score,a)#HR
  return(c(summary(cc)$coefficients[1,c(1,3)],bb$concordance,bb$std.err))
}

modelMetaStaticsMULTI<-function(model_file,input_file){
  model_score<-scale(modelScore(model_file,input_file))
  clinical_model<-clinical[[which(input_file==names(clinical))]]
  a<-c()
  if(input_file!="Patch et al."){
    nonNA_po<-match(rownames(model_score),rownames(clinical_model))
    a<-cbind(clinical_model$days_to_death[nonNA_po],
             clinical_model$vital_status[nonNA_po],model_score,
             clinical_model$tumorstage,clinical_model$grade)
    a[,2]<-ifelse(a[,2]=="living",0,1)
  }else{
    nonNA_po<-match(rownames(model_score),rownames(clinical_model))
    a<-cbind(clinical_model$donor_survival_time[nonNA_po],
             clinical_model$donor_vital_status[nonNA_po],model_score,
             clinical_model$tumour_stage,clinical_model$tumour_grade)
    a[,2]<-ifelse(a[,2]=="alive",0,1)
  }
  rownames(a)<-names(model_score)
  colnames(a)<-c("time","status","model_score","stage","grade")
  if(idx.tumorstage[i]&(!idx.tumorgrade[i])){
    a<-a[,-5]
  }else if((!idx.tumorstage[i])&(!idx.tumorgrade[i])){
    a<-a[,-c(4:5)]
  }else if((!idx.tumorstage[i])&(idx.tumorgrade[i])){
    a<-a[,-4]
  }
  if(length(unique(a[,])))
  a<-as.data.frame(apply(a,2,as.numeric))
  bb<-survConcordance(Surv(time,status)~model_score,a)#一致性
  cc<-coxph(Surv(time,status)~.,a)#HR
  return(c(summary(cc)$coefficients[1,c(1,3)],bb$concordance,bb$std.err))
}
#####modelMetaStatics-auc#####
modelMetaStatics.auc<-function(model_score,input_file,time){
  #model_score<-scale(modelScore(model_file,input_file))
  clinical_model<-clinical[[which(input_file==names(clinical))]]
  a<-c()
  if(input_file!="Patch et al."){
    nonNA_po<-match(names(model_score),rownames(clinical_model))
    a<-cbind(clinical_model$days_to_death[nonNA_po],
             clinical_model$vital_status[nonNA_po],model_score)
    a[,2]<-ifelse(a[,2]=="living",0,1)
  }else{
    nonNA_po<-match(names(model_score),rownames(clinical_model))
    a<-cbind(clinical_model$donor_survival_time[nonNA_po],
             clinical_model$donor_vital_status[nonNA_po],model_score)
    a[,2]<-ifelse(a[,2]=="alive",0,1)
  }
  rownames(a)<-names(model_score)
  colnames(a)<-c("time","status","model_score")
  phe<-as.data.frame(apply(a,2,as.numeric))
  
  a<-timeROC(T=phe$time,
             delta=phe$status,marker=phe$model_score,
             cause=1,weighting = "marginal",
             times=time,iid = TRUE)
  b<-c()
  if(length(time)==2){
    b<-c(a$AUC[1],a$inference$vect_sd_1[1],a$AUC[2],a$inference$vect_sd_1[2])
    names(b)<-c("aucOf3","se.aucof3","aucOf5","se.aucof5")
  }else{
    b<-c(a$AUC[2],a$inference$vect_sd_1[2])
    names(b)<-c("aucOf3","se.aucof3")
  }
  
  return(b)
}