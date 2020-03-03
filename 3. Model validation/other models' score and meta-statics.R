source('/code/lower2symmetric.R')
source('/codes/signatureSimility.R')
load("model_coefs1.Rdata")
load("datasets_scale.Rdata")
globalExp<-data.frame(GN=rownames(datasets_scale[[1]]),datasets_scale[[1]])
for(i in 2:length(datasets_scale)){
  globalExp<-merge(globalExp,data.frame(GN=rownames(datasets_scale[[i]]),
                                        datasets_scale[[i]]),
                   by.x="GN",by.y="GN",all=T)
}
rownames(globalExp)<-globalExp$GN
globalExp<-globalExp[,-1]
ModelGeneSim<-c()
ModelGeneSimP<-c()
for(j in 1:13){
  for(i in (j+1):14){
    sa<-as.character(names(model.coefs1[[i]]))
    sb<-as.character(names(model.coefs1[[j]]))
    ModelGeneSim<-c(ModelGeneSim,signatureSimility(sa,sb,globalExp))
    ramsim<-c()
    for(m in 1:100){
      set.seed(m);ram.a<-sample(rownames(globalExp),length(sa))
      set.seed(m+1);ram.b<-sample(rownames(globalExp),length(sb))
      ramsim[m]<-signatureSimility(ram.a,ram.b,globalExp)
    }
    z_score<-as.numeric(scale(ramsim))
    SD<-sd(ramsim)#
    MEAN<-mean(ramsim)
    real.scale<-as.numeric((ModelGeneSim[length(ModelGeneSim)]-MEAN)/SD)
    ModelGeneSimP<-c(ModelGeneSimP,2*pnorm(abs(real.scale),lower.tail=FALSE))
  }
}
ModelGeneSim<-lower2symmetric(14,ModelGeneSim)
ModelGeneSimP<-lower2symmetric(16,ModelGeneSimP)

input_files<-names(datasets_nonscale)[trainPo]
model_files<-names(model.coefs1)
predictive_ability_Train<-list()
for(i in 1:length(input_files)){
  predictive_ability_Train.i<-c()
  for(j in 1:length(model.coefs1)){
    predictive_ability_Train.i<-rbind(predictive_ability_Train.i,
                                      modelMetaStatics(model_files[j],input_files[i]))
    cat(i,j,"\n")
  }
  rownames(predictive_ability_Train.i)<-names(model.coefs1)
  colnames(predictive_ability_Train.i) <- c("ln(HR)", "se(ln(HR))",
                                            "C-index","se(C-index)")
  predictive_ability_Train[[i]]<-predictive_ability_Train.i
}
names(predictive_ability_Train)<-names(datasets_scale)[trainPo]

predictive_ability_Train_excel<-c()
for(i in 1:length(predictive_ability_Train)){
  predictive_ability_Train_excel<-rbind(predictive_ability_Train_excel,predictive_ability_Train[[i]])
}
predictive_ability_Train_excel<-Digits4(predictive_ability_Train_excel)
colnames(predictive_ability_Train_excel)<-c("lnHR","se lnHR","C-index","se C-index")
changeNamePo<-seq(from=1,to=1+14*(length(predictive_ability_Train)-1),by=14)
rownames(predictive_ability_Train_excel)[changeNamePo]<-paste(names(datasets_nonscale)[trainPo],
                                                              rownames(predictive_ability_Train_excel)[changeNamePo],
                                                              sep="%")
                                                              
matrix_meta_Train<-c()
for(i in 1:nrow(predictive_ability_Train[[1]])){
  matrix_premeta<-c()
  for(j in 1:length(predictive_ability_Train)){
    matrix_premeta<-rbind(matrix_premeta,predictive_ability_Train[[j]][i,])
  }
  lnHR<-meta.summaries(matrix_premeta[,1],matrix_premeta[,2],method="fixed",logscale=T)
  if(lnHR$het[3]>0.05){
    summary_HR<-c(lnHR$summary,lnHR$se.summary)
  }else{
    lnHR<-meta.summaries(matrix_premeta[,1],matrix_premeta[,2],method="random",logscale=T)
    summary_HR<-c(lnHR$summary,lnHR$se.summary)
  }
  C<-meta.summaries(matrix_premeta[,3],matrix_premeta[,4],method="fixed",logscale=F)
  if(C$het[3]>0.05){
    summary_C<-c(C$summary,C$se.summary)
  }else{
    C<-meta.summaries(matrix_premeta[,3],matrix_premeta[,4],method="random",logscale=F)
    summary_C<-c(C$summary,C$se.summary)
  }
  matrix_meta_Train<-rbind(matrix_meta_Train,c(summary_HR,summary_C))
}
rownames(matrix_meta_Train)<-rownames(predictive_ability_Train[[1]])
colnames(matrix_meta_Train)<-c("lnHR","se(lnHR)","C-index","se(C-index)")

datasets_scale_train<-datasets_scale[trainPo]
datasets_nonscale_train<-datasets_nonscale[trainPo]
clinicalTrain<-clinical[trainPo]
clinicalMatrixTrain<-getSurvivalMatrix(clinicalTrain,"global","Patch et al.")
clinicalMatrixTrain[which("Alive"==clinicalMatrixTrain[,2]),2]<-0
clinicalMatrixTrain[which("Dead"==clinicalMatrixTrain[,2]),2]<-1
ForteenPerfoTrain<-c()
input_files<-names(datasets_nonscale)[trainPo]
model_files<-names(model.coefs1)
for(i in 1:length(model.coefs1)){
  modelGlobalScore.i<-c()
  for(j in 1:length(input_files)){
    modelGlobalScore.i<-c(modelGlobalScore.i,
                          modelScore(model_files[i],input_files[j]))
  }
  modelGlobalScore.i<-scale(modelGlobalScore.i)
  if(!is.null(modelGlobalScore.i)){
    globalDataTrain<-data.frame(time=as.numeric(clinicalMatrixTrain[,1]),
                                status=as.numeric(clinicalMatrixTrain[,2]),
                                model_score=modelGlobalScore.i)
    bb<-survConcordance(Surv(time,status)~model_score,globalDataTrain)#一致性
    cc<-coxph(Surv(time,status)~model_score,globalDataTrain)#HR
    ForteenPerfoTrain<-rbind(ForteenPerfoTrain,c(summary(cc)$coefficients[1,c(1,3)],bb$concordance,bb$std.err))
    rownames(ForteenPerfoTrain)[nrow(ForteenPerfoTrain)]<-names(model.coefs1)[i]
  }
}
colnames(ForteenPerfoTrain)<-c("lnHR","selnHR","C-index","se C-index")
input_files<-names(datasets_nonscale)[-trainPo]
model_files<-names(model.coefs1)
predictive_ability_Test<-list()
for(i in 1:length(input_files)){
  predictive_ability_Test.i<-c()
  for(j in 1:length(model.coefs1)){
    predictive_ability_Test.i<-rbind(predictive_ability_Test.i,
                                     modelMetaStatics(model_files[j],input_files[i]))
    cat(i,j,"\n")
  }
  rownames(predictive_ability_Test.i)<-names(model.coefs1)
  colnames(predictive_ability_Test.i) <- c("ln(HR)", "se(ln(HR))",
                                           "C-index","se(C-index)")
  predictive_ability_Test[[i]]<-predictive_ability_Test.i
}
names(predictive_ability_Test)<-names(datasets_scale)[-trainPo]
save(predictive_ability_Test,file="predictive_ability_Test.Rdata")
predictive_ability_Test_excel<-c()
for(i in 1:length(predictive_ability_Test)){
  predictive_ability_Test_excel<-rbind(predictive_ability_Test_excel,predictive_ability_Test[[i]])
}
predictive_ability_Test_excel<-Digits4(predictive_ability_Test_excel)
colnames(predictive_ability_Test_excel)<-c("lnHR","se lnHR","C-index","se C-index")
changeNamePo<-seq(from=1,to=1+14*(length(predictive_ability_Test)-1),by=14)
rownames(predictive_ability_Test_excel)[changeNamePo]<-paste(names(datasets_nonscale)[-trainPo],
                                                             rownames(predictive_ability_Test_excel)[changeNamePo],
                                                             sep="%")
                                                             
matrix_meta_Test<-c()
for(i in 1:nrow(predictive_ability_Test[[1]])){
  matrix_premeta<-c()
  for(j in 1:length(predictive_ability_Test)){
    matrix_premeta<-rbind(matrix_premeta,predictive_ability_Test[[j]][i,])
  }
  lnHR<-meta.summaries(matrix_premeta[,1],matrix_premeta[,2],method="fixed",logscale=T)
  if(lnHR$het[3]>0.05){
    summary_HR<-c(lnHR$summary,lnHR$se.summary)
  }else{
    lnHR<-meta.summaries(matrix_premeta[,1],matrix_premeta[,2],method="random",logscale=T)
    summary_HR<-c(lnHR$summary,lnHR$se.summary)
  }
  C<-meta.summaries(matrix_premeta[,3],matrix_premeta[,4],method="fixed",logscale=F)
  if(C$het[3]>0.05){
    summary_C<-c(C$summary,C$se.summary)
  }else{
    C<-meta.summaries(matrix_premeta[,3],matrix_premeta[,4],method="random",logscale=F)
    summary_C<-c(C$summary,C$se.summary)
  }
  matrix_meta_Test<-rbind(matrix_meta_Test,c(summary_HR,summary_C))
}
rownames(matrix_meta_Test)<-rownames(predictive_ability_Test[[1]])
colnames(matrix_meta_Test)<-c("lnHR","se(lnHR)","C-index","se(C-index)")

datasets_scale_test<-datasets_scale[-trainPo]
datasets_nonscale_test<-datasets_nonscale[-trainPo]
clinicalTest<-clinical[-trainPo]
clinicalMatrixTest<-getSurvivalMatrix(clinicalTest,"global","Patch et al.")
clinicalMatrixTest[which("Alive"==clinicalMatrixTest[,2]),2]<-0
clinicalMatrixTest[which("Dead"==clinicalMatrixTest[,2]),2]<-1
ForteenPerfoTest<-c()
input_files<-names(datasets_nonscale)[-trainPo]
model_files<-names(model.coefs1)
for(i in 1:length(model.coefs1)){
  modelGlobalScore.i<-c()
  for(j in 1:length(input_files)){
    modelGlobalScore.i<-c(modelGlobalScore.i,
                          modelScore(model_files[i],input_files[j]))
    name<-c(name,names(modelGlobalScore.i))
  }
  name<-names(modelGlobalScore.i)
  modelGlobalScore.i<-scale(modelGlobalScore.i)
  names(modelGlobalScore.i)<-name
  if(!is.null(modelGlobalScore.i)){
    globalDataTest<-data.frame(time=as.numeric(clinicalMatrixTest[name,1]),
                               status=as.numeric(clinicalMatrixTest[name,2]),
                               model_score=modelGlobalScore.i)
    bb<-survConcordance(Surv(time,status)~model_score,globalDataTest)#一致性
    cc<-coxph(Surv(time,status)~model_score,globalDataTest)#HR
    ForteenPerfoTest<-rbind(ForteenPerfoTest,c(summary(cc)$coefficients[1,c(1,3)],bb$concordance,bb$std.err))
    rownames(ForteenPerfoTest)[nrow(ForteenPerfoTest)]<-names(model.coefs1)[i]
  }
}
colnames(ForteenPerfoTest)<-c("lnHR","selnHR","C-index","se C-index")