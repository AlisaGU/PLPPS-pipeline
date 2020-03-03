load("clinical.Rdata")
load("ssgsea_out.Rdata")
load("ssgsea_out_test.Rdata")
load("ssgsea_out_all.Rdata")
load("ssgsea.out.independ.Rdata")
load("datasets_scale.Rdata")
clinical.chemotherapy<-list()
#####train####
{
  ctrainpo<-9
  preclinical<-t(apply(as.matrix(clinical[[ctrainpo]][,31],ncol=1), 1,function(x){unlist(strsplit(x,split = "///"))}))
  rownames(preclinical)<-rownames(clinical[[ctrainpo]])
  colnames(preclinical)<-unlist(lapply(strsplit(preclinical[1,],": "), function(x){ifelse(length(x)==2,x[1],x[2])}))
  preclinical<-as.data.frame(t(apply(preclinical,1,function(x){unlist(lapply(strsplit(x,": "), function(x){ifelse(length(x)==2,x[2],x[3])}))})))
  clinical.chemotherapy[["GSE32062"]]<-data.frame(score=structure(ssgsea.out[[4]][3,],
                                                                  .Names=colnames(ssgsea.out[[4]])),
  )}

ctrainpo<-13
a<-data.frame(score=structure(ssgsea.out[[7]][3,],
                              .Names=colnames(ssgsea.out[[7]])),
              fufa.status=clinical[[13]]$recurrence_status,
              fufa.days=clinical[[13]]$days_to_tumor_recurrence,
              os.status=clinical[[13]]$vital_status,
              os.days=clinical[[13]]$days_to_death,
              chemo.agent=NA,
              chemo.response=clinical[[13]]$primary_therapy_outcome_success)
clinical.chemotherapy[["Dressman"]]<-a
#####test####
ctestpo<-8
preclinical<-t(apply(as.matrix(clinical[[ctestpo]][,31],ncol=1), 1,function(x){unlist(strsplit(x,split = "///"))}))
rownames(preclinical)<-rownames(clinical[[ctestpo]])
colnames(preclinical)<-unlist(lapply(strsplit(preclinical[1,],": "), function(x){ifelse(length(x)==2,x[1],x[2])}))
preclinical<-as.data.frame(t(apply(preclinical,1,function(x){unlist(lapply(strsplit(x,": "), function(x){ifelse(length(x)==2,x[2],x[3])}))})))
write.table(preclinical,"clinical.chemo.8.txt",quote = F,sep="\t")
clinical.chemo.8<-read.table("clinical.chemo.8.txt",header = T,sep="\t",as.is = T)
a<-data.frame(score=structure(ssgsea.out.test[[5]][3,],
                              .Names=colnames(ssgsea.out.test[[5]])),
              fufa.status=clinical[[8]][,17],
              fufa.days=clinical[[8]][,16],
              os.status=clinical[[8]]$vital_status,
              os.days=clinical[[8]]$days_to_death,
              chemo.agent=clinical.chemo.8[,23],
              chemo.response=clinical.chemo.8[,22])
a<-a[-which(a$chemo.response=="unknown"),]
clinical.chemotherapy[["GSE30161"]]<-a

{
  ctestpo<-10
  preclinical<-t(apply(as.matrix(clinical[[ctestpo]][,31],ncol=1), 1,function(x){unlist(strsplit(x,split = "///"))}))
  rownames(preclinical)<-rownames(clinical[[ctestpo]])
  colnames(preclinical)<-unlist(lapply(strsplit(preclinical[1,],": "), function(x){ifelse(length(x)==2,x[1],x[2])}))
  preclinical<-as.data.frame(t(apply(preclinical,1,function(x){unlist(lapply(strsplit(x,": "), function(x){ifelse(length(x)==2,x[2],x[3])}))})))}

ctest_independpo<-14
a<-data.frame(score=structure(ssgsea.out.independ[[1]][3,],
                              .Names=colnames(ssgsea.out.independ[[1]])),
              fufa.status=clinical[[14]]$RFS,
              fufa.days=clinical[[14]]$RFS.time,
              os.status=clinical[[14]]$vital_status,
              os.days=clinical[[14]]$days_to_death,
              chemo.agent=clinical[[14]]$additional_pharmaceutical_therapy,
              chemo.response=clinical[[14]]$primary_therapy_outcome_success)
clinical.chemotherapy[["TCGA"]]<-a[-which(a$chemo.response==""),]


ctest_independpo<-16
a<-data.frame(score=structure(ssgsea.out.independ[[2]][3,],
                              .Names=colnames(ssgsea.out.independ[[2]])),
              fufa.status=clinical[[16]]$disease_status_last_followup,
              fufa.days=clinical[[16]]$donor_relapse_interval,
              os.status=clinical[[16]]$donor_vital_status,
              os.days=clinical[[16]]$donor_survival_time,
              chemo.agent=clinical[[16]]$specimen_donor_treatment_type_other,
              chemo.response=clinical[[16]]$Response)
clinical.chemotherapy[["patch"]]<-a

library(ggplot2)
source('G:/haodapeng-data/deleteScatter/code/chemoLabel.R')



chemo<-rbind(rbind(rbind(chemo1,chemo2),chemo3),chemo4)
chemo<-chemo[order(chemo$score),]

label<-rep(1:43,each=10)
a<-c()
for(i in 1:10){
  x<-chemo[(43*(i-1)+1):(43*i),7]
  a<-c(a,length(c(which(x=="CR (completely response)"),
                  which(x=="Sensitive"),
                  which(x=="Complete Remission/Response"),
                  which(x=="completeresponse")))/length(x))
}
Data<-data.frame(phase=1:10,
                 CR.propo=a)

ggplot(Data,aes(x=phase,y=CR.propo))+
  geom_point()+
  xlim(c(0,10))+
  labs(y="Complete response (%)",x="Phase")+
  geom_smooth(method = lm,se=F)+
  theme_bw()+
  theme(
    axis.text.x = element_text(color="black",size = 12),
    axis.text.y = element_text(color="black",size = 12),
    axis.title.x = element_text(color="black",size = 16),
    legend.text = element_text(color="black",size = 12),
    legend.title =element_text(color="black",size = 12),
    panel.border = element_blank(),panel.grid=element_blank(),
    axis.line = element_line(colour = "black")
  )

fit<-lm(CR.propo~phase,Data)
summary(fit)
