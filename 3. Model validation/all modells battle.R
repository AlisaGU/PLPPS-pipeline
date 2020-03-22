
HR_C_Train<-read.csv("HR_CTrain.csv",header=T,row.names = 1,as.is = T)
load("predictive_ability_Train.Rdata")
pdf("14 VS ssgsea seperate train.pdf",width = 7,height = 7)
for(i in 1:length(predictive_ability_Train)){
  ForModel<-predictive_ability_Train[[i]]
  ssgsea<-HR_C_Train[i,c(2,3,6,7)]
  colnames(ForModel)<-names(ssgsea)<-c("lnHR","selnHR","C_index","se_C_index")
  Data<-rbind(ForModel,ssgsea)
  rownames(Data)<-c(rownames(Data)[1:14],"ssgsea")
  
  Data<-as.data.frame(Data[order(Data[,1],decreasing = T),])
  meta.i<-metagen(lnHR,selnHR,sm="HR",data = Data,
                  studlab = rownames(Data))
  forest(meta.i, leftcols = c("studlab"),leftlabs = c("Study Element"),
         rightcols = c("effect.ci"),fontsize = 13,just.studlab = "center",
         overall = F,col.square = "#363886",lwd=2,
         xlab=names(predictive_ability_Train)[i],
         weight.study = "same",col.inside="black",ref = 1)
  
  Data<-as.data.frame(Data[order(Data[,3],decreasing = T),])
  meta.i<-metagen(C_index,se_C_index,sm="C-index",data = Data,
                  studlab = rownames(Data))
  forest(meta.i, leftcols = c("studlab"),leftlabs = c("Study Element"),
         rightcols = c("effect.ci"),fontsize = 13,just.studlab = "center",
         overall = F,col.square = "#363886",lwd=2,
         xlab=names(predictive_ability_Train)[i],
         weight.study = "same",col.inside="black",ref = 0.5)
}
dev.off()

ssgsea_colla_Train<-as.numeric(unlist(read.csv("ssgsea_colla_Train.csv",
                                               header=T,row.names = 1,as.is = T)))
matrix_meta_Train<-read.csv("matrix_meta_Train.csv",header=T,row.names = 1,as.is = T)#ѵ???ϲ?
colnames(matrix_meta_Train)<-names(ssgsea_colla_Train)<-c("lnHR","selnHR",
                                                          "C_index","se_C_index")

Data<-rbind(matrix_meta_Train,ssgsea_colla_Train)
rownames(Data)<-c(rownames(matrix_meta_Train),"-FIPI")


Data<-as.data.frame(Data[order(Data[,3],decreasing = T),])
meta.i<-metagen(C_index,se_C_index,sm="C-index",data = Data,
                studlab = rownames(Data))
forest(meta.i, leftcols = c("studlab"),leftlabs = c("Study Element"),
       rightcols = c("effect.ci"),fontsize = 13,just.studlab = "center",
       overall = F,col.square = "#363886",lwd=2,xlim = c(0.45,0.7),plotwidth="3.5cm",
       weight.study = "same",col.inside="black",ref = 0.5,squaresize = 0.5)
       


Data<-as.data.frame(Data[order(Data[,1],decreasing = T),])
meta.i<-metagen(lnHR,selnHR,sm="lnHR",data = Data,
                studlab = rownames(Data))
forest(meta.i, leftcols = c("studlab"),leftlabs = c("Study Element"),
       rightcols = c("effect.ci"),fontsize = 13,just.studlab = "center",
       overall = F,col.square = "#363886",lwd=2,plotwidth="3.5cm",squaresize = 0.5,
       weight.study = "same",col.inside="black",ref = 0,xlim = c(-0.2,0.6))
       
       
colnames(global_Value_se)<-colnames(ForteenPerfoTrain)<-c("lnHR","selnHR",
                                                          "C_index","se_C_index")
Data<-rbind(ForteenPerfoTrain,global_Value_se[1,])
rownames(Data)<-c(rownames(ForteenPerfoTrain),"-FIPI")
pdf("14 VS ssgsea globalTrainC-index.pdf",width=7,height = 7)
Data<-as.data.frame(Data[order(Data[,3],decreasing = T),])
meta.i<-metagen(C_index,se_C_index,sm="C-index",data = Data,
                studlab = rownames(Data))
forest(meta.i, leftcols = c("studlab"),leftlabs = c("Study Element"),
       rightcols = c("effect.ci"),fontsize = 13,just.studlab = "center",
       overall = F,col.square = "#363886",lwd=2,xlim = c(0.4,0.7),
       weight.study = "same",col.inside="black",ref = 0.5)
dev.off()


Data<-as.data.frame(Data[order(Data[,1],decreasing = T),])
meta.i<-metagen(lnHR,selnHR,sm="HR",data = Data,
                studlab = rownames(Data))
forest(meta.i, leftcols = c("studlab"),leftlabs = c("Study Element"),
       rightcols = c("effect.ci"),fontsize = 13,just.studlab = "center",
       overall = F,col.square = "#363886",lwd=2,xlim = c(0.9,7),
       weight.study = "same",col.inside="black",ref = 1)
       

HR_C_Test<-read.csv("HR_C_Test.csv",header=T,row.names = 1,as.is = T)
load("predictive_ability_Test.Rdata")

for(i in 1:length(predictive_ability_Test)){
  ForModel<-predictive_ability_Test[[i]]
  ssgsea<-HR_C_Test[i,c(2,3,6,7)]
  colnames(ForModel)<-names(ssgsea)<-c("lnHR","selnHR","C_index","se_C_index")
  Data<-rbind(ForModel,ssgsea)
  rownames(Data)<-c(rownames(Data)[1:14],"ssgsea")
  
  Data<-as.data.frame(Data[order(Data[,1],decreasing = T),])
  meta.i<-metagen(lnHR,selnHR,sm="HR",data = Data,
                  studlab = rownames(Data))
  forest(meta.i, leftcols = c("studlab"),leftlabs = c("Study Element"),
         rightcols = c("effect.ci"),fontsize = 13,just.studlab = "center",
         overall = F,col.square = "#363886",lwd=2,
         xlab=names(predictive_ability_Train)[i],
         weight.study = "same",col.inside="black",ref = 1)
  
  Data<-as.data.frame(Data[order(Data[,3],decreasing = T),])
  meta.i<-metagen(C_index,se_C_index,sm="C-index",data = Data,
                  studlab = rownames(Data))
  forest(meta.i, leftcols = c("studlab"),leftlabs = c("Study Element"),
         rightcols = c("effect.ci"),fontsize = 13,just.studlab = "center",
         overall = F,col.square = "#363886",lwd=2,
         xlab=names(predictive_ability_Train)[i],
         weight.study = "same",col.inside="black",ref = 0.5)
}

ssgsea_colla_Test<-as.numeric(unlist(read.csv("ssgsea_colla_Test.csv",
                                              header=T,row.names = 1,as.is = T)))
matrix_meta_Test<-read.csv("matrix_meta_Test.csv",header=T,row.names = 1,as.is = T)#ѵ???ϲ?
colnames(matrix_meta_Test)<-names(ssgsea_colla_Test)<-c("lnHR","selnHR",
                                                        "C_index","se_C_index")

Data<-rbind(matrix_meta_Test,ssgsea_colla_Test)
rownames(Data)<-c(rownames(ForteenPerfoTest),"ssgsea")


Data<-as.data.frame(Data[order(Data[,3],decreasing = T),])
meta.i<-metagen(C_index,se_C_index,sm="C-index",data = Data,
                studlab = rownames(Data))
forest(meta.i, leftcols = c("studlab"),leftlabs = c("Study Element"),
       rightcols = c("effect.ci"),fontsize = 13,just.studlab = "center",
       overall = F,col.square = "#363886",lwd=2,xlim = c(0.4,0.7),
       weight.study = "same",col.inside="black",ref = 0.5)
       


Data<-as.data.frame(Data[order(Data[,1],decreasing = T),])
meta.i<-metagen(lnHR,selnHR,sm="HR",data = Data,
                studlab = rownames(Data))
forest(meta.i, leftcols = c("studlab"),leftlabs = c("Study Element"),
       rightcols = c("effect.ci"),fontsize = 13,just.studlab = "center",
       overall = F,col.square = "#363886",lwd=2,
       weight.study = "same",col.inside="black",ref = 1)
       
colnames(global_Value_se)<-colnames(ForteenPerfoTest)<-c("lnHR","selnHR",
                                                         "C_index","se_C_index")
Data<-rbind(ForteenPerfoTest,global_Value_se[2,])
rownames(Data)<-c(rownames(ForteenPerfoTest),"ssgsea")

Data<-as.data.frame(Data[order(Data[,3],decreasing = T),])
meta.i<-metagen(C_index,se_C_index,sm="C-index",data = Data,
                studlab = rownames(Data))
forest(meta.i, leftcols = c("studlab"),leftlabs = c("Study Element"),
       rightcols = c("effect.ci"),fontsize = 13,just.studlab = "center",
       overall = F,col.square = "#363886",lwd=2,xlim = c(0.4,0.7),
       weight.study = "same",col.inside="black",ref = 0.5)
       


Data<-as.data.frame(Data[order(Data[,1],decreasing = T),])
meta.i<-metagen(lnHR,selnHR,sm="HR",data = Data,
                studlab = rownames(Data))
forest(meta.i, leftcols = c("studlab"),leftlabs = c("Study Element"),
       rightcols = c("effect.ci"),fontsize = 13,just.studlab = "center",
       overall = F,col.square = "#363886",lwd=2,xlim = c(0.9,3),
       weight.study = "same",col.inside="black",ref = 1)
       
       
load("model_coefs1.Rdata")
jaccardMatrix_model<-matrix(nrow=14,ncol=14)
colnames(jaccardMatrix_model)<-names(model.coefs1)
rownames(jaccardMatrix_model)<-names(model.coefs1)
for(i in 1:14){
  for(j in 1:14){
    intersect.ij<-length(intersect(as.character(names(model.coefs1[[i]])),
                                   as.character(names(model.coefs1[[j]]))))
    union.ij<-length(union(names(model.coefs1[[i]]),names(model.coefs1[[j]])))
    jaccard<-intersect.ij/union.ij
    jaccardMatrix_model[i,j]<-jaccard
  }
}
diag(jaccardMatrix_model)<-NA


data<-data.frame(class=rep(rownames(jaccardMatrix_model),each=13),
                 value=unlist(jaccardMatrix_model)[!is.na(unlist(jaccardMatrix_model))])
ggplot(data,aes(class,value))+
  geom_boxplot(fill="#EFC000",outlier.shape = NA)+
  geom_jitter(shape=21,fill="#EFC000",size=1.7,width=0.3,height=0)+
  theme_bw()+
  labs(y="Jaccard Index")+
  scale_y_continuous(breaks=c(0.00,0.01,0.02),labels=c("0.00","0.01","0.02"))+
  theme(panel.border = element_blank(),panel.grid=element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1,color = "black",
                                   size=14),
        axis.text.y=element_text(color="black",size=14),
        plot.title = element_text(colour = "black", face = "bold", 
                                  size = 14, vjust = 1,hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y=element_text(color="black",size=16),
        legend.position ="none")
        


gradeStatus<-c(0,rep(1,5),0,rep(1,8),2)
gradeStatus[14]<-2
result<-HGClinicalStatus(clinical,gradeStatus)
idx.tumorstage<-result[,1]
idx.tumorgrade<-result[,2]
idx.age<-result[,3]
stageName<-c(rep("tumorstage",15),"tumour_stage")
gradeName<-c(rep("grade",15),"tumour_grade")
ageName<-c(rep("age_at_initial_pathologic_diagnosis",15),
           "donor_age_at_diagnosis")
timeName<-c(rep("days_to_death",15),"donor_survival_time")
eventName<-c(rep("vital_status",15),"donor_vital_status")
selfEvalResult<-list()
for(i in 1:length(clinical)){
  score.i<-scale(as.numeric(ssgsea.out.all[[i]][3,]))
  phe<-clinical[[i]]
  
  sCol<-which(colnames(phe)==stageName[i])[1]
  gCol<-which(colnames(phe)==gradeName[i])[1]
  aCol<-which(colnames(phe)==ageName[i])[1]
  timeCol<-which(colnames(phe)==timeName[i])[1]
  eventCol<-which(colnames(phe)==eventName[i])[1]
  
  result<-selfEval(score.i,phe,names(clinical)[i],
                   idx.tumorstage[i],idx.tumorgrade[i],
                   idx.age[i],sCol,gCol,aCol,timeCol,eventCol)
  selfEvalResult[[names(clinical)[i]]]<-result
}
uni_ssgsea_colla<-unicoxForestData(selfEvalResult)[[1]]
multi_ssgsea_colla<-MulticoxForestData(selfEvalResult)[[1]]

#####others#####
#####uni####
input_files<-names(datasets_nonscale)
model_files<-names(model.coefs1)
uni_predictive_ability<-list()
for(i in 1:length(input_files)){
  predictive_ability.i<-c()
  for(j in 1:length(model.coefs1)){
    predictive_ability.i<-rbind(predictive_ability.i,
                                modelMetaStatics(model_files[j],input_files[i]))
    cat(i,j,"\n")
  }
  rownames(predictive_ability.i)<-names(model.coefs1)
  colnames(predictive_ability.i) <- c("ln(HR)", "se(ln(HR))",
                                      "C-index","se(C-index)")
  uni_predictive_ability[[i]]<-predictive_ability.i
}
names(uni_predictive_ability)<-names(datasets_scale)

uni_matrix_meta<-c()
for(i in 1:nrow(uni_predictive_ability[[1]])){
  matrix_premeta<-c()
  for(j in 1:length(uni_predictive_ability)){
    matrix_premeta<-rbind(matrix_premeta,uni_predictive_ability[[j]][i,])
  }
  lnHR<-metagen(matrix_premeta[,1],matrix_premeta[,2],sm="lnHR")
  if(lnHR$pval.Q>0.05){
    summary_HR<-c(lnHR$TE.fixed,lnHR$seTE.fixed,lnHR$pval.fixed)
  }else{
    summary_HR<-c(lnHR$TE.random,lnHR$seTE.random,lnHR$pval.random)
  }
  C<-metagen(matrix_premeta[,3],matrix_premeta[,4],sm="C-index",null.effect=0.5 )
  if(C$pval.Q>0.05){
    summary_C<-c(C$TE.fixed,C$seTE.fixed,C$pval.fixed)
  }else{
    summary_C<-c(C$TE.random,C$seTE.random,C$pval.random)
  }
  uni_matrix_meta<-rbind(uni_matrix_meta,c(summary_HR,summary_C))
}
rownames(uni_matrix_meta)<-rownames(uni_predictive_ability[[1]])
colnames(uni_matrix_meta)<-c("lnHR","se(lnHR)","p.HR",
                             "C-index","se(C-index)","p.C")
#####multi#####
input_files<-names(datasets_nonscale)
model_files<-names(model.coefs1)
multi_predictive_ability<-list()
for(i in 1:length(input_files)){
  predictive_ability.i<-c()
  for(j in 1:length(model.coefs1)){
    predictive_ability.i<-rbind(predictive_ability.i,
                                modelMetaStaticsMULTI(model_files[j],input_files[i]))
    cat(i,j,"\n")
  }
  rownames(predictive_ability.i)<-names(model.coefs1)
  colnames(predictive_ability.i) <- c("ln(HR)", "se(ln(HR))",
                                      "C-index","se(C-index)")
  multi_predictive_ability[[i]]<-predictive_ability.i
}
names(multi_predictive_ability)<-names(datasets_scale)

Multi_matrix_meta<-c()
for(i in 1:nrow(multi_predictive_ability[[1]])){
  matrix_premeta<-c()
  for(j in 1:length(multi_predictive_ability)){
    matrix_premeta<-rbind(matrix_premeta,multi_predictive_ability[[j]][i,])
  }
  lnHR<-metagen(matrix_premeta[,1],matrix_premeta[,2],sm="lnHR")
  if(lnHR$pval.Q>0.05){
    summary_HR<-c(lnHR$TE.fixed,lnHR$seTE.fixed,lnHR$pval.fixed)
  }else{
    summary_HR<-c(lnHR$TE.random,lnHR$seTE.random,lnHR$pval.random)
  }
  C<-metagen(matrix_premeta[,3],matrix_premeta[,4],sm="C-index",null.effect=0.5)
  if(C$pval.Q>0.05){
    summary_C<-c(C$TE.fixed,C$seTE.fixed,C$pval.fixed)
  }else{
    summary_C<-c(C$TE.random,C$seTE.random,C$pval.random)
  }
  Multi_matrix_meta<-rbind(Multi_matrix_meta,c(summary_HR,summary_C))
}
rownames(Multi_matrix_meta)<-rownames(multi_predictive_ability[[1]])
colnames(Multi_matrix_meta)<-c("lnHR","se(lnHR)","p.HR",
                               "C-index","se(C-index)","p.C")
#####1orest#####
#####uni####
Data<-rbind(uni_matrix_meta,uni_ssgsea_colla)
rownames(Data)<-c(rownames(uni_matrix_meta),"PLPPS")
rownames(Data)[which("ours116"==rownames(Data))]<-"lu14"
colnames(Data)<-c("lnHR","se(lnHR)","p.HR",
                  "C-index","se(C-index)","p.C")
Data<-as.data.frame(Data[order(Data[,1],decreasing = T),])
write.csv(Data,"uniHR_C_p.csv",quote=F)

{
  Data<-read.csv("uniHR_C_p.csv",header=T,row.names = 1,as.is=T)
  Data<-as.data.frame(Data[order(Data[,1],decreasing = F),])
  Data$p.HR<--log10(Data$p.HR)
  Data$model<-factor(rownames(Data),levels=rownames(Data))
  ggplot(Data,aes(x=model,y=p.HR))+
    geom_bar(position=position_dodge(), stat="identity")+
    xlab("")+coord_flip()+
    scale_y_continuous(expand=c(0,0))+
    ylab(expression(paste(-log[10],"(p)")))+
    theme_bw()+
    labs(fill="")+
    theme(panel.border = element_blank(),panel.grid=element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text.x = element_text(size=15,color = "black"),
          axis.text.y=element_text(size=15,color="black"),
          plot.title = element_text(colour = "black", face = "bold", 
                                    size = 14, vjust = 1,hjust = 0.5),
          axis.title.x = element_text(size=15,color = "black"),
          axis.title.y = element_text(size=15,color = "black"),
          legend.text = element_text(size = 15),
          legend.position="top")

{
  Data<-read.csv("uniHR_C_p.csv",header=T,row.names = 1,as.is=T)
  Data<-as.data.frame(Data[order(Data[,1],decreasing = F),])
  Data$p.C<--log10(Data$p.C)
  Data$model<-factor(rownames(Data),levels=rownames(Data))
  ggplot(Data,aes(x=model,y=p.C))+
    geom_bar(position=position_dodge(), stat="identity")+
    xlab("")+coord_flip()+
    scale_y_continuous(expand=c(0,0))+
    ylab(expression(paste(-log[10],"(p)")))+
    theme_bw()+
    labs(fill="")+
    theme(panel.border = element_blank(),panel.grid=element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text.x = element_text(size=15,color = "black"),
          axis.text.y=element_text(size=15,color="black"),
          plot.title = element_text(colour = "black", face = "bold", 
                                    size = 14, vjust = 1,hjust = 0.5),
          axis.title.x = element_text(size=15,color = "black"),
          axis.title.y = element_text(size=15,color = "black"),
          legend.text = element_text(size = 15),
          legend.position="top")}

uniHRname<-rownames(Data)
meta.i<-metagen(Data[,1],Data[,2],sm="HR",data = Data,
                studlab = rownames(Data))
forest(meta.i, leftcols = c("studlab"),leftlabs = c("Study Element"),
       rightcols = c("effect.ci"),fontsize = 13,just.studlab = "center",
       at=c(0.9,1,1.1,1.2,1.3,1.4,1.5,1.6),
       overall = F,col.square = "#363886",lwd=2,squaresize = 0.4,xlim=c(0.9,1.6),
       weight.study = "same",col.inside="black",ref = 1,plotwidth="3.5cm",
       hetlab ="",print.I2=F,print.tau2=F,print.pval.Q=F)


meta.i<-metagen(Data[,4],Data[,5],sm="C-index",data = Data,
                studlab = rownames(Data))
forest(meta.i, leftcols = c("studlab"),leftlabs = c("Study Element"),
       rightcols = c("effect.ci"),fontsize = 13,just.studlab = "center",
       overall = F,col.square = "#363886",lwd=2,xlim = c(0.4,0.7),
       at=c(0.4,0.5,0.6,0.7),plotwidth="3.5cm",squaresize = 0.4,
       weight.study = "same",col.inside="black",ref = 0.5,
       hetlab ="",print.I2=F,print.tau2=F,print.pval.Q=F)

Data<-Data[c(15,3,1,6,5,2,7),]
Data<-as.data.frame(Data[order(Data[,3],decreasing = T),])
Data$model<-factor(rownames(Data),levels =rownames(Data) )
ggplot(Data, aes(x=model, y=C.index)) + 
  geom_bar(stat ="identity",width = 0.8,fill="#BE0005") +
  geom_errorbar(aes(ymin=C.index-1.96*se.C.index, ymax=C.index+1.96*se.C.index), width=.1) +
  labs(y="C-index")+
  scale_y_continuous(expand = c(0,0),breaks = seq(0,0.65,0.05))+
  coord_cartesian(ylim = c(0.4, 0.7))+
  theme_bw()+
  theme(panel.border = element_blank(),panel.grid=element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1,color = "black",
                                   size=14),
        axis.text.y=element_text(color="black",size=14),
        plot.title = element_text(colour = "black", face = "bold", 
                                  size = 14, vjust = 1,hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y=element_text(color="black",size=16),
        legend.position ="none")

#####multi#####
Data<-rbind(Multi_matrix_meta,multi_ssgsea_colla)
rownames(Data)<-c(rownames(Multi_matrix_meta),"PLPPS")
rownames(Data)[which("ours116"==rownames(Data))]<-"lu14"
colnames(Data)<-c("lnHR","se(lnHR)","p.HR",
                  "C-index","se(C-index)","p.C")

Data<-as.data.frame(Data[order(Data[,1],decreasing = T),])

{
  Data<-as.data.frame(Data[order(Data[,1],decreasing = F),])
  Data$p.HR<--log10(Data$p.HR)
  Data$model<-factor(rownames(Data),levels=rownames(Data))
  ggplot(Data,aes(x=model,y=p.HR))+
    geom_bar(position=position_dodge(), stat="identity")+
    xlab("")+coord_flip()+
    scale_y_continuous(expand=c(0,0))+
    ylab(expression(paste(-log[10],"(p)")))+
    theme_bw()+
    labs(fill="")+
    theme(panel.border = element_blank(),panel.grid=element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text.x = element_text(size=15,color = "black"),
          axis.text.y=element_text(size=15,color="black"),
          plot.title = element_text(colour = "black", face = "bold", 
                                    size = 14, vjust = 1,hjust = 0.5),
          axis.title.x = element_text(size=15,color = "black"),
          axis.title.y = element_text(size=15,color = "black"),
          legend.text = element_text(size = 15),
          legend.position="top")}

MultiHRname<-rownames(Data)
meta.i<-metagen(lnHR,selnHR,sm="HR",data = Data,
                studlab = rownames(Data))
forest(meta.i, leftcols = c("studlab"),leftlabs = c("Study Element"),
       rightcols = c("effect.ci"),fontsize = 13,just.studlab = "center",
       at=c(0.9,1,1.1,1.2,1.3,1.4,1.6),
       overall = F,col.square = "#363886",lwd=2,squaresize = 0.4,xlim=c(0.9,1.6),
       weight.study = "same",col.inside="black",ref = 1,plotwidth="3.5cm")

Data<-as.data.frame(Data[order(Data[,3],decreasing = T),])
MultiCname<-rownames(Data)
meta.i<-metagen(C.index,se.C.index,sm="C-index",data = Data,
                studlab = rownames(Data))
forest(meta.i, leftcols = c("studlab"),leftlabs = c("Study Element"),
       rightcols = c("effect.ci"),fontsize = 13,just.studlab = "center",
       overall = F,col.square = "#363886",lwd=2,xlim = c(0.4,0.7),
       at=c(0.4,0.5,0.6,0.7),plotwidth="3.5cm",squaresize = 0.4,
       weight.study = "same",col.inside="black",ref = 0.5)
#####1heatmap#####
#####uni#####
HR_C<-unicoxForestData(selfEvalResult)[[2]]
k<-HR_C$lnHR
o<-c()
for(i in 1:length(uni_predictive_ability)){
  o<-rbind(o,uni_predictive_ability[[i]][,1])
}
Data<-cbind(k,o)
rownames(Data)<-rownames(HR_C)
colnames(Data)[c(1,15)]<-c("PLPPS","lu14")
bk<-c(min(Data),seq(-0.5,-0.01,length.out = 49),
      0,seq(0.01,0.5,length.out = 49),max(Data))
color<-c("#7A7ABC",colorRampPalette(c("#7A7ABC","white"))(48),"white",
         colorRampPalette(c("white","firebrick3"))(49),"firebrick3")
Data<-t(Data)
Data<-Data[uniHRname,]
pheatmap(Data,cluster_cols = F,cluster_rows = F,
         color=color,breaks=bk,border_color="white",
         cellwidth = 15,cellheight = 15)


k<-HR_C$C_index
o<-c()
for(i in 1:length(uni_predictive_ability)){
  o<-rbind(o,uni_predictive_ability[[i]][,3])
}
Data<-cbind(k,o)
rownames(Data)<-rownames(HR_C)
colnames(Data)[c(1,15)]<-c("PLPPS","lu14")
bk<-c(seq(0.3,0.49,length.out=20),0.5,seq(0.51,0.7,length.out=20),max(Data))
color<-c(colorRampPalette(c("#7A7ABC","white"))(19),"white",
         colorRampPalette(c("white","firebrick3"))(20),"firebrick3")
Data<-t(Data)
Data<-Data[uniHRname,]
pheatmap(Data,cluster_cols = F,cluster_rows = F,
         color=color,breaks=bk,border_color="white",
         cellwidth = 15,cellheight = 15)

#####multi####
HR_C<-MulticoxForestData(selfEvalResult)[[2]]
k<-HR_C$lnHR
o<-c()
for(i in 1:length(multi_predictive_ability)){
  o<-rbind(o,multi_predictive_ability[[i]][,1])
}
Data<-cbind(k,o)
rownames(Data)<-rownames(HR_C)
colnames(Data)[c(1,15)]<-c("-FIPI","lu14")
Data<-t(Data)
Data<-Data[MultiHRname,]
bk<-c(min(Data),seq(-0.5,-0.01,length.out=49),0,seq(0.01,0.5,length.out=49),max(Data))
color<-c("#7A7ABC",colorRampPalette(c("#7A7ABC","white"))(48),"white",
         colorRampPalette(c("white","firebrick3"))(49),"firebrick3")
pheatmap(Data,cluster_cols = F,cluster_rows = F,
         color=color,breaks=bk,border_color="white",
         cellwidth = 15,cellheight = 15)


k<-HR_C$C_index
o<-c()
for(i in 1:length(multi_predictive_ability)){
  o<-rbind(o,multi_predictive_ability[[i]][,3])
}
Data<-cbind(k,o)
rownames(Data)<-rownames(HR_C)
colnames(Data)[c(1,15)]<-c("-FIPI","lu14")
bk<-c(seq(0.35,0.49,length.out=15),0.5,seq(0.51,0.7,length.out=20),max(Data))
color<-c(colorRampPalette(c("#7A7ABC","white"))(15),"white",
         colorRampPalette(c("white","firebrick3"))(20),"firebrick3")
pheatmap(Data,cluster_cols = F,cluster_rows = F,
         color=color,breaks=bk,border_color="white",
         cellwidth = 15,cellheight = 15)
       