#####persetting#####
source('/code/signatureScore.R')
source('/code/getTotalGene.R')
source('/code/ROC_CI_global.R')
source('/code/globalEval.R')
source('/code/HGClinicalStatus.R')
source('/code/plotSurvivalRoc.R')
source('/code/sstselfEvalUMContinue.R')
source('/code/selfEval.R')
source('/code/sstselfEvalUM.R')
source('/code/unicoxForestData.R')
source('/code/MulticoxForestData.R')
source('/code/sstselfEval.R')
source('/code/groupQ.R')
source('/code/createTable.R')
source('/code/ROC_CI.R')
source('/code/combineCoxResult.R')
source('/code/getSurvivalMatrix.R')
source('/code/plotSurvivalCurve.R')
source('/code/getSurvivalMatrix1.R')
source('/code/group.R')
source('/code/Digits4.R')
source('/code/getGlobalExp.R')
source('/code/getGlobalScore.R')
source('/code/getHeatmapExp.R')
source('/code/ROC3_5.R')
source('/code/ROC3_5_CI.R')

library(rmeta)
library(ggplot2)
library(survminer)
library(survival)
library(survivalROC)
library(GEOquery)
library(Hmisc)
library(meta)
library(clusterProfiler)
library(pheatmap)
library(timeROC)
load("model_method.RData")
load("model_coefs1.Rdata")
load("datasets_nonscale.Rdata")
load("datasets_scale.Rdata")
load("ssgsea_out.Rdata")
load("ssgsea_out_test.Rdata")
load("ssgsea_out_all.Rdata")
load("ssgsea.out.independ.Rdata")
load("clinical.Rdata")
trainPo<-which(sapply(datasets_scale,ncol)>100)[-8]
testPo<-which(sapply(datasets_scale,ncol)<100)[-8]
independPo<-c(14,16)
riskModelGene[22]<-"SPART"
score<-structure(as.numeric(unlist(sapply(ssgsea.out,function(x){x[3,]}))),
                 .Names=as.character(unlist(sapply(ssgsea.out,colnames))))
medianScore<-median(score)
Q1score<--4786.0
Q3score<- -906.4
gene<-c(protectModelGene,riskModelGene)
score.test<-structure(as.numeric(unlist(sapply(ssgsea.out.test,function(x){x[3,]}))),
                      .Names=as.character(unlist(sapply(ssgsea.out.test,colnames))))
score.all<-structure(as.numeric(unlist(sapply(ssgsea.out.all,function(x){x[3,]}))),
                     .Names=as.character(unlist(sapply(ssgsea.out.all,colnames))))
score.test.independ<-structure(as.numeric(unlist(sapply(ssgsea.out.independ,function(x){x[3,]}))),
                               .Names=as.character(unlist(sapply(ssgsea.out.independ,colnames))))

clinicalTrain<-clinical[trainPo]
clinicalTest<-clinical[testPo]
clinicalIndepend<-clinical[independPo]
clinicalMatrixTrain<-getSurvivalMatrix(clinicalTrain,"global","Patch et al.")
clinicalMatrixTest<-getSurvivalMatrix(clinicalTest,"global","Patch et al.")

x1<-clinical[[independPo[1]]]
phe1<-cbind(x1$days_to_death,x1$vital_status)
x2<-clinical[[independPo[[2]]]]
phe2<-cbind(x2$donor_survival_time,x2$donor_vital_status)
phe2[,2]<-ifelse(phe2[,2]==2,1,0)
phe<-rbind(phe1,phe2)
rownames(phe)<-c(rownames(x1),rownames(x2))
colnames(phe)<-c("time","status")
phe[,2]<-ifelse(phe[,2]=="0","Alive","Dead")
clinicalMatrixIndepend<-phe

clinicalMatrixTestAll<-rbind(clinicalMatrixTest,clinicalMatrixIndepend)
clinicalMatrixAll<-rbind(clinicalMatrixTrain,
                         clinicalMatrixTest,
                         clinicalMatrixIndepend)

allROC<-list()

model.coef1.cox<-model.coefs1[c(1,2,3,5,6,11)]
model_files<-names(model.coef1.cox)
uniHRname<-c("PLPPS","yoshihara12","tcga11",
             "crijns09","yoshihara10","bonome08_263",
             "lu14","mok09","kernagis12",
             "bonome08_572","denkert09","sabatier11",
             "kang12","konstantinopoulos10","hernandez10")
#####training#####
datasets_scale_train<-datasets_scale[trainPo]
datasets_nonscale_train<-datasets_nonscale[trainPo]
input_files<-names(datasets_nonscale)[trainPo]

model7score.train<-list()
for(i in 1:length(model.coef1)){
  modelGlobalScore.i<-list()
  for(j in 1:length(input_files)){
    modelGlobalScore.i[[input_files[j]]]<-modelScore(model_files[i],input_files[j])
  }
  model7score.train[[names(model.coef1)[i]]]<-modelGlobalScore.i
  
}
model7score.train[["ssgsea"]]<-score

ppp<-c()
for(i in 1:length(model7score.train)){
  ppp<-rbind(ppp,ROC3_5_CI(scale(model7score.train[[i]]),clinicalMatrixTrain,1,2))
  
}
rownames(ppp)<-names(model7score.train)
colnames(ppp)<-c("3year","5year")
allROC[["ssgsea_cox"]][["train"]]<-ppp


data<-data.frame(Model=rep(rownames(ppp),2),
                 AUC=c(ppp[,1],ppp[,2]),
                 Group=c(rep("3 years",7),
                         rep("5 years",7)))
input_files<-names(datasets_nonscale)
model_files<-c(names(model.coefs1),"ssgsea")

ROC<-list()
for(i in 1:length(input_files)){
  ROC.i<-c()
  time<-c()
  if(i==3 |i==11){
    time<-c(365.25*3)
  }else{
    time<-c(365.25*3,365.25*5)
  }
  for(j in 1:length(model_files)){
    modelscore.ij<-c()
    if(j==15){
      modelscore.ij<-structure(ssgsea.out.all[[i]][3,],
                               names=colnames(ssgsea.out.all[[i]]))
    }else{
      modelscore.ij<-modelScore(model_files[j],input_files[i])
    }
    
    ROC.i<-rbind(ROC.i,modelMetaStatics.auc(modelscore.ij,input_files[i],time))
    cat(i,j,"\n")
  }
  rownames(ROC.i)<-model_files
  ROC[[i]]<-ROC.i
}
names(ROC)<-input_files
y3_matrix_meta<-c()
for(i in 1:nrow(ROC[[1]])){
  matrix_premeta<-c()
  for(j in 1:length(ROC)){
    matrix_premeta<-rbind(matrix_premeta,ROC[[j]][i,1:2])
  }
  AUC<-metagen(matrix_premeta[,1],matrix_premeta[,2],sm="AUC",null.effect=0.5)
  if(AUC$pval.Q>0.05){
    summary_AUC<-c(AUC$TE.fixed,AUC$seTE.fixed,AUC$pval.fixed)
  }else{
    summary_AUC<-c(AUC$TE.random,AUC$seTE.random,AUC$pval.random)
  }
  y3_matrix_meta<-rbind(y3_matrix_meta,summary_AUC)
}

y5_matrix_meta<-c()
for(i in 1:nrow(ROC[[1]])){
  matrix_premeta<-c()
  for(j in c(1:length(ROC))[-c(3,11)]){
    matrix_premeta<-rbind(matrix_premeta,ROC[[j]][i,3:4])
  }
  AUC<-metagen(matrix_premeta[,1],matrix_premeta[,2],sm="AUC",null.effect=0.5)
  if(AUC$pval.Q>0.05){
    summary_AUC<-c(AUC$TE.fixed,AUC$seTE.fixed,AUC$pval.fixed)
  }else{
    summary_AUC<-c(AUC$TE.random,AUC$seTE.random,AUC$pval.random)
  }
  y5_matrix_meta<-rbind(y5_matrix_meta,summary_AUC)
}


Multi_matrix_meta<-cbind(y3_matrix_meta,y5_matrix_meta)

rownames(Multi_matrix_meta)<-rownames(ROC[[1]])
colnames(Multi_matrix_meta)<-c("aucOf3","se.aucof3","p.aucof3",
                               "aucOf5","se.aucof5","p.aucof5")

Data<-Multi_matrix_meta
rownames(Data)[15]<-"PLPPS"
rownames(Data)[14]<-"lu14"
Data<-as.data.frame(Data[match(uniHRname,rownames(Data)),])

write.csv(Data[,1:3],"AUC_p.csv",quote=F)
{
  Data<-read.csv("AUC_p.csv",header=T,row.names = 1,as.is=T)
  Data<-Data[15:1,]
  Data$p.aucof3<--log10(Data$p.aucof3)
  Data$model<-factor(rownames(Data),levels=rownames(Data))
  ggplot(Data,aes(x=model,y=p.aucof3))+
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
  ggsave("AUC_p.pdf",width = 7,height=7)}

Data$model<-factor(rownames(Data),levels =rownames(Data))
meta.i<-metagen(aucOf3,se.aucof3,sm="AUC",data = Data,
                studlab = rownames(Data),null.effect = 0.5)
pdf("allmodels'forest_metaROC.pdf",width=7,height = 7)
forest(meta.i, leftcols = c("studlab"),leftlabs = c("Study Element"),
       rightcols = c("effect.ci"),fontsize = 13,just.studlab = "center",
       overall = F,col.square = "#363886",lwd=2,xlim = c(0.4,0.7),
       at=c(0.4,0.5,0.6,0.7),plotwidth="3.5cm",squaresize = 0.4,
       weight.study = "same",col.inside="black",ref = 0.5,
       hetlab ="",print.I2=F,print.tau2=F,print.pval.Q=F)
dev.off()

heatmapData<-c()
for(i in 1:length(ROC)){
  heatmapData<-cbind(heatmapData,ROC[[i]][,1])
}
colnames(heatmapData)<-names(ROC)
rownames(heatmapData)[14]<-"lu14"
rownames(heatmapData)[15]<-"PLPPS"
heatmapData<-heatmapData[match(rownames(Data),rownames(heatmapData)),]
bk<-c(seq(0.3,0.49,length.out=19),0.5,seq(0.51,0.81,length.out=30))
color<-c(colorRampPalette(c("#7A7ABC","white"))(19),"white",
         colorRampPalette(c("white","firebrick3"))(30),"firebrick3")
pdf("allmodels_meta_ROC-heatmap.pdf",width=6,height=5)
pheatmap(heatmapData,cluster_cols = F,cluster_rows = F,
         color=color,breaks=bk,border_color="white",
         cellwidth = 15,cellheight = 15)
dev.off()



Data<-Multi_matrix_meta[c(1,2,3,5,6,7,15),]
Data<-as.data.frame(Data[order(Data[,1],decreasing = T),])
Data$model<-factor(rownames(Data),levels =rownames(Data) )
ggplot(Data, aes(x=model, y=aucOf3)) + 
  geom_bar(stat ="identity",width = 0.8,fill="#BE0005") +
  geom_errorbar(aes(ymin=aucOf3-1.96*se.aucof3, ymax=aucOf3+1.96*se.aucof3), width=.1) +
  labs(y="Survival ROC")+
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
ggsave("coxmodels_meta_ROC.pdf",width = 7,height=7)