source('/code/getSurvivalMatrix.R')
source('/code/getSurvivalMatrix1.R')
source('/code/mergeexp.R')
library(survival)
library(survminer)
library(ggplot2)
library(DMwR)
trainPo<-which(sapply(datasets_scale,ncol)>100)[-8]
testPo<-which(sapply(datasets_scale,ncol)<100)[-8]
independPo<-c(14,16)
oriexp<-datasets_nonscale
modelgene<-c(protectModelGene,riskModelGene)
#####training#####
clinicalTrain<-clinical[trainPo]
clinicalMatrix<-getSurvivalMatrix(clinicalTrain,"global","Patch et al.")
clinicalMatrix[,2]<-ifelse(clinicalMatrix[,2]=="Dead",1,0)
list<-oriexp[trainPo]
trainexp<-matrix(NA,nrow=68,ncol=1190)
rownames(trainexp)<-modelgene
colnames(trainexp)<-unlist(sapply(list, colnames))
for(i in 1:length(list)){
  rowpo.list<-match(rownames(trainexp),rownames(list[[i]]))
  colpo<-match(colnames(list[[i]]),colnames(trainexp))
  trainexp[,colpo]<-list[[i]][rowpo.list,]
}
trainexp<-knnImputation(as.data.frame(trainexp),scale = F)
trainexp1<-as.data.frame(t(trainexp))
multicox<-coxph(Surv(as.numeric(clinicalMatrix[,1]),
                     as.numeric(clinicalMatrix[,2]))~.,data = trainexp1)
multicox.sum<-summary(multicox)
gene68_coxmodel<-structure(multicox.sum$coefficients[,1],
                           names=rownames(multicox.sum$coefficients))
names(gene68_coxmodel)[c(12,25)]<-c("HLA-E","HLA-F")

#score
Data<-trainexp
score<-structure(apply(Data,2,function(x){return(sum((x*gene68_coxmodel)))}),
                 names=colnames(Data))
MEDIAN<-median(score)
group<-ifelse(score>=MEDIAN,"HighRisk","LowRisk")
survivalMatrix<-as.data.frame(cbind(clinicalMatrix,score,group))
fit<-survfit(Surv(as.numeric(as.character(survivalMatrix$time)),
                  as.numeric(as.character(survivalMatrix$status)))~group,
             data=survivalMatrix)
diff<-survdiff(Surv(as.numeric(as.character(survivalMatrix$time)),
                    as.numeric(as.character(survivalMatrix$status)))~group,
               data=survivalMatrix)
ggsurvplot(fit, data = survivalMatrix, title="Meta-training (n = 1190)",
           risk.table = TRUE,risk.table.col = "strata",
           pval=T,
           risk.table.height = 0.2,xlab="Time(Days)",legend.title="",
           legend.labs=c("High Risk","Low Risk"),risk.table.title="No. at risk",
           palette = c("#E7201C","#0A499E"),
           ggtheme = theme_bw()+theme(panel.border = element_blank(),
                                      panel.grid=element_blank(),
                                      axis.line = element_line(colour = "black"),
                                      axis.text.x = element_text(size=15,color = "black"),
                                      axis.text.y=element_text(size=15,color="black"),
                                      plot.title = element_text(colour = "black", face = "bold", 
                                                                size = 14, vjust = 1,hjust = 0.5),
                                      axis.title.x = element_text(size=15,color = "black"),
                                      axis.title.y = element_text(size=15,color = "black"),
                                      legend.text = element_text(size = 15),
                                      legend.title = element_text(size = 15)))
#####est1#####
clinicalTest<-clinical[testPo]
clinicalMatrix<-getSurvivalMatrix(clinicalTest,"global","Patch et al.")
clinicalMatrix[,2]<-ifelse(clinicalMatrix[,2]=="Dead",1,0)

list<-oriexp[testPo]
test1exp<-matrix(NA,nrow=68,ncol=438)
rownames(test1exp)<-modelgene
colnames(test1exp)<-unlist(sapply(list, colnames))
for(i in 1:length(list)){
  rowpo.list<-match(rownames(test1exp),rownames(list[[i]]))
  colpo<-match(colnames(list[[i]]),colnames(test1exp))
  test1exp[,colpo]<-list[[i]][rowpo.list,]
}
test1exp<-knnImputation(as.data.frame(test1exp),scale = F)
score.test<-apply(test1exp,2,function(x){return(sum((x*gene68_coxmodel)))})
group.test<-ifelse(score.test>=MEDIAN,"HighRisk","LowRisk")
survivalMatrix<-as.data.frame(cbind(clinicalMatrix,score.test,group.test))
fit<-survfit(Surv(as.numeric(as.character(survivalMatrix$time)),
                  as.numeric(as.character(survivalMatrix$status)))~group.test,
             data=survivalMatrix)
diff<-survdiff(Surv(as.numeric(as.character(survivalMatrix$time)),
                    as.numeric(as.character(survivalMatrix$status)))~group.test,
               data=survivalMatrix)
ggsurvplot(fit, data = survivalMatrix, title="Meta-testing I(n = 438)",
           risk.table = TRUE,risk.table.col = "strata",
           #pval =  format(1-pchisq(diff$chisq,1),digits = 3),
           pval=T,
           risk.table.height = 0.2,xlab="Time(Days)",legend.title="",
           legend.labs=c("High Risk","Low Risk"),risk.table.title="No. at risk",
           palette = c("#E7201C","#0A499E"),
           ggtheme = theme_bw()+theme(panel.border = element_blank(),panel.grid=element_blank(),
                                      axis.line = element_line(colour = "black"),
                                      axis.text.x = element_text(size=15,color = "black"),
                                      axis.text.y=element_text(size=15,color="black"),
                                      plot.title = element_text(colour = "black", face = "bold", 
                                                                size = 14, vjust = 1,hjust = 0.5),
                                      axis.title.x = element_text(size=15,color = "black"),
                                      axis.title.y = element_text(size=15,color = "black"),
                                      legend.text = element_text(size = 15),
                                      legend.title = element_text(size = 15)))
#####test2#####
clinicalIndepend<-clinical[independPo]#TCGA、patch
clinicalMatrix11<-getSurvivalMatrix1(clinicalIndepend[1],"global","Patch et al.")
clinicalMatrix22<-getSurvivalMatrix(clinicalIndepend[2],"global","Patch et al.")
clinicalMatrix<-rbind(clinicalMatrix11,clinicalMatrix22)
clinicalMatrix[,2]<-ifelse(clinicalMatrix[,2]=="Dead",1,0)


list<-oriexp[independPo]
test2exp<-matrix(NA,nrow=68,ncol=354)
rownames(test2exp)<-modelgene
colnames(test2exp)<-unlist(sapply(list, colnames))
for(i in 1:length(list)){
  rowpo.list<-match(rownames(test2exp),rownames(list[[i]]))
  colpo<-match(colnames(list[[i]]),colnames(test2exp))
  test2exp[,colpo]<-as.matrix(list[[i]][rowpo.list,])
}
score.test1<-apply(test2exp,2,function(x){return(sum((x*gene68_coxmodel)))})
group.test1<-ifelse(score.test1>=MEDIAN,"HighRisk","LowRisk")
survivalMatrix<-as.data.frame(cbind(clinicalMatrix,score.test1,group.test1))
fit<-survfit(Surv(as.numeric(as.character(survivalMatrix$time)),
                  as.numeric(as.character(survivalMatrix$status)))~group.test1,
             data=survivalMatrix)
diff<-survdiff(Surv(as.numeric(as.character(survivalMatrix$time)),
                    as.numeric(as.character(survivalMatrix$status)))~group.test1,
               data=survivalMatrix)
ggsurvplot(fit, data = survivalMatrix, title="Meta-testing II(n = 354)",
           risk.table = TRUE,risk.table.col = "strata",
           #pval =  format(1-pchisq(diff$chisq,1),digits = 3),
           pval=T,
           risk.table.height = 0.2,xlab="Time(Days)",legend.title="",
           legend.labs=c("High Risk","Low Risk"),risk.table.title="No. at risk",
           palette = c("#E7201C","#0A499E"),
           ggtheme = theme_bw()+theme(panel.border = element_blank(),panel.grid=element_blank(),
                                      axis.line = element_line(colour = "black"),
                                      axis.text.x = element_text(size=15,color = "black"),
                                      axis.text.y=element_text(size=15,color="black"),
                                      plot.title = element_text(colour = "black", face = "bold", 
                                                                size = 14, vjust = 1,hjust = 0.5),
                                      axis.title.x = element_text(size=15,color = "black"),
                                      axis.title.y = element_text(size=15,color = "black"),
                                      legend.text = element_text(size = 15),
                                      legend.title = element_text(size = 15)))
#画图
data<-data.frame(pvalue=-log10(c(0.000364,0.00359,0.2,0.01)),
                 Method=c(rep("PLPPS",2),rep("cox",2)),
                 dataset=factor(rep(c("Meta-testing I","Meta-testing II"),2),
                                levels = c("Meta-testing I","Meta-testing II")))
ggplot(data,aes(x=dataset,y=pvalue,fill=Method))+
  geom_bar(position=position_dodge(), stat="identity")+
  xlab("")+
  # scale_x_continuous(breaks=c(1,2,3),
  #                    labels=c("Meta-testing I","Meta-testing II"))+
  scale_y_continuous(expand=c(0,0),limits = c(0,4))+
  ylab(expression(paste(-log[10],"(p)")))+
  scale_fill_manual(values =c("#E16253","#7CA7DE"))+
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
