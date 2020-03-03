#####TCGA-seq独立数据集处理######
setwd("G:\\haodapeng-data\\deleteScatter\\think1\\test")
HiSeqV2<-read.table("HiSeqV2",header=T,row.names = 1,as.is=T,sep="\t")
OV_clinicalMatrix<-read.table("OV_clinicalMatrix",header=T,as.is=T,sep="\t")
OV_clinicalMatrix$sampleID<-gsub("-",".",OV_clinicalMatrix$sampleID,fixed=T)
stage<-(OV_clinicalMatrix$clinical_stage=="Stage IIIA")|(OV_clinicalMatrix$clinical_stage=="Stage IIIB")|(OV_clinicalMatrix$clinical_stage=="Stage IIIC")|(OV_clinicalMatrix$clinical_stage=="Stage IV")
grade<-(OV_clinicalMatrix$neoplasm_histologic_grade=="G2")|OV_clinicalMatrix$neoplasm_histologic_grade=="G3"
OV_clinicalMatrix1<-OV_clinicalMatrix[(!is.na(OV_clinicalMatrix$X_EVENT))&(!is.na(OV_clinicalMatrix$OS.time))&stage&grade,]
o<-intersect(OV_clinicalMatrix1$sampleID,colnames(HiSeqV2))
exp<-HiSeqV2[,o]
OV_clinicalMatrix2<-OV_clinicalMatrix1[match(o,OV_clinicalMatrix1$sampleID),]
colnames(OV_clinicalMatrix2)[which(colnames(OV_clinicalMatrix2)=="X_EVENT")]<-"vital_status"
colnames(OV_clinicalMatrix2)[which(colnames(OV_clinicalMatrix2)=="OS.time")]<-"days_to_death"
colnames(OV_clinicalMatrix2)[which(colnames(OV_clinicalMatrix2)=="neoplasm_histologic_grade")]<-"grade"
colnames(OV_clinicalMatrix2)[which(colnames(OV_clinicalMatrix2)=="clinical_stage")]<-"tumorstage"
TCGA.seq<-list(exp=exp,phe=OV_clinicalMatrix2)
save(TCGA.seq,file="TCGA_seq.Rdata")
save(TCGA.seq,file="G:\\haodapeng-data\\deleteScatter\\think1\\TCGA_seq.Rdata")

setwd("G:\\haodapeng-data\\deleteScatter\\think1")
load("clinical.Rdata")
load("datasets_scale.Rdata")
load("datasets_nonscale.Rdata")
clinical$TCGA<-OV_clinicalMatrix2
datasets_nonscale$TCGA<-TCGA.seq[[1]]
a<-t(scale(t(TCGA.seq[[1]])))
datasets_scale$TCGA<-a
save(clinical,file="clinical.Rdata")
save(datasets_nonscale,file="datasets_nonscale.Rdata")
save(datasets_scale,file="datasets_scale.Rdata")
#####数据集详细信息表#####
options(digits=3)
library(GEOquery)
setwd("G:\\haodapeng-data\\deleteScatter\\think1")
source('G:/haodapeng-data/deleteScatter/code/getSurvivalMatrix1.R')
source('G:/haodapeng-data/deleteScatter/code/getSurvivalMatrix.R')
source('G:/haodapeng-data/deleteScatter/code/stage_grade_age.R')
load("datasets_scale.Rdata")
load("clinical.Rdata")
dataset<-names(datasets_scale)
use<-rep(NA,16)
trainPo<-which(sapply(clinical,nrow)>100)[-8]
testPo<-which(sapply(clinical,nrow)<100)[-8]
independPo<-c(14,16)
use[trainPo]<-"Training"
use[testPo]<-"Test"
use[independPo]<-"Independent Test"

MedianSurvivalandCensor<-c()
for(i in 1:length(clinical)){
  if(i!=14){
    clinical.i<-getSurvivalMatrix(clinical[[i]],type="seperate",
                                  exceptSet="Patch et al.",DataName = names(clinical)[i])
  }else{
    clinical.i<-getSurvivalMatrix1(clinical[[i]],type="seperate",
                                  exceptSet="Patch et al.",DataName = names(clinical)[i])
  }
  MedianSurvivalandCensor<-rbind(MedianSurvivalandCensor,
                                 c(format(median(as.numeric(clinical.i[,1])),digits=4),
                                   format(length(which(clinical.i[,2]=="Dead"))*100/nrow(clinical.i),digits=2)))
}
platform<-c("Illumina HumanRef-8 v2","Operon human v3","Affymetrix  U133a",
            "Agilent G4112F", "Affymetrix U133 Plus 2.0","Affymetrix U133 Plus 2.0",
            "Affymetrix  U133a", "Affymetrix U133 Plus 2.0","Agilent G4112F",
            "Agilent G4112F","ABI Human Genome Survey Microarray Version 2",
            "Affymetrix U133 Plus 2.0","Affymetrix  U133a","RNA seq",
            "Agilent G4110B","RNA seq")
Accession<-c("22348002","19192944","19294737",
             "20300634","19962670","22101765",
             "18593951","22348014","22241791",
             "22241791","22497737","18698038",
             "17290060","21720365","24368280","26017449")
totalSample<-c()#所有样本数
lateSerousSample<-c()#晚期浆液性卵巢癌样本数
overallType<-rep("OS",16)#生存时间类型
#sex<-c()#男性的数目和比例
STAGE<-c()#stage>=3的样本数目和比例
GRADE<-c()#grade>=3的样本数目和比例
AGE<-c()#年龄的均值和95%置信区间
gradeStatus<-c(0,rep(1,5),0,rep(1,8),2)
gradeStatus[14]<-2
stageColname<-c(rep("tumorstage",15),"tumour_stage")
gradeColname<-c(rep("grade",15),"tumour_grade")
ageColname<-c(rep("age_at_initial_pathologic_diagnosis",15),
              "donor_age_at_diagnosis")
clinical[[14]]$tumorstage<-ifelse(grepl("III",clinical[[14]]$tumorstage),3,4)
clinical[[14]]$grade<-ifelse(grepl("G2",clinical[[14]]$grade),2,3)
for(i in 1:16){
  phe<-clinical[[i]]
  totalSample<-c(totalSample,nrow(phe))
  result<-stage_grade_age(phe,gradeStatus = gradeStatus[i],
                          stageColname = stageColname[i],
                          gradeColname = gradeColname[i],
                          ageColname = ageColname[i])
  lateSerousSample<-c(lateSerousSample,result[1])
  STAGE<-c(STAGE,result[2])
  GRADE<-c(GRADE,result[3])
  AGE<-c(AGE,result[4])
}
detailsOfDataset<-data.frame(dataset,use,platform,totalSample,
                             lateSerousSample,overallType,STAGE,
                             GRADE,AGE,MedianSurvivalandCensor,Accession)
colnames(detailsOfDataset)<-c("Dataset","Use","Platform","Number of OV patients",
                              "Number of HGSC patients","Overall Type",
                              "Stage >T3","Grade >=G3","Age(95% CI,year)",
                              "Median Survival(day)","Censored(%)","Ref(PMID)")
write.csv(detailsOfDataset,"detailsOfDataset.csv",row.names = F)
save(clinical,file="clinical.Rdata")
#####用cox获取生存相关基因####
setwd("G:\\haodapeng-data\\deleteScatter\\think1")
source('G:/haodapeng-data/deleteScatter/code/HGClinicalStatus.R')
source('G:/haodapeng-data/deleteScatter/code/coxAnalysis.R')
library(GEOquery)
library(survival) 
load("clinical.Rdata")
load("HGSC datasets with OS.rda")
gradeStatus<-c(0,rep(1,5),0,rep(1,8),2)
gradeStatus[14]<-2
#获取各个数据集HG样本的stage,grade,age状态
result<-HGClinicalStatus(clinical,gradeStatus)
idx.tumorstage<-result[,1]
idx.tumorgrade<-result[,2]
idx.age<-result[,3]
expNa<-c(rep(0,13),1,1)
# load("gene_cox.Rdata")
# load("gene_meta.Rdata")
# load("datasets_scale.Rdata")
#获取cox生存相关基因
gene.cox<-list()
gene.meta<-list()
datasets_scale<-list()
datasets_nonscale<-list()
clinical<-list()
for(i in 1:15){
  MyGSE <- HDP.OVsets[[i]];phe <- pData(MyGSE);Mat <- exprs(MyGSE)
  if(gradeStatus[i]==0){
    phe<-phe[which(phe$histological_type=="ser" & phe$sample_type=="tumor" & phe$summarystage=="late"),]
  }else if(gradeStatus[i]==1){
    phe<-phe[which(phe$histological_type=="ser" & phe$sample_type=="tumor" & phe$summarystage=="late" & phe$grade>=2),]
  }
  clinical[[names(HDP.OVsets)[i]]]<-phe
  Mat <- Mat[,rownames(phe)];datasets_nonscale[[names(HDP.OVsets)[i]]]<-Mat
  Mat <- t(scale(t(Mat)));datasets_scale[[names(HDP.OVsets)[i]]]<-Mat
  result.cox<-coxAnalysis(Mat,phe,idx.tumorstage[i],
                          idx.tumorgrade[i],idx.age[i],
                          expNa[i])
  gene.cox[[names(HDP.OVsets)[i]]]<-result.cox[[1]]
  gene.meta[[names(HDP.OVsets)[i]]]<-result.cox[[2]]
  gene.cox[[14]]<-result.cox[[1]]
  gene.meta[[14]]<-result.cox[[2]]
}

i=16;MyGSE <- HDP.OVsets[[i]];phe <- pData(MyGSE);
clinical[[names(HDP.OVsets)[i]]]<-phe;
Mat <- exprs(MyGSE);
Mat<-Mat[-which(apply(Mat,1,function(x){all(x==x[1])})),]#删掉的是基因
datasets_nonscale[[names(HDP.OVsets)[16]]]<-Mat
Mat <- t(scale(t(Mat)))
datasets_scale[[names(HDP.OVsets)[16]]]<-Mat
a<-c()
d<-c()
for (j in 1:nrow(Mat)) {
  gene.j <- coxph(Surv(donor_survival_time,donor_vital_status == "deceased")
                  ~ Mat[j,]+tumour_stage+tumour_grade+
                    donor_age_at_diagnosis, phe)
  b<-summary(gene.j)
  a<-rbind(a,c(b$coefficients[1,5],b$conf.int[1,c(1,3,4)]))
  d<-rbind(d,b$coefficients[1,1:3])
}
rownames(a)<-rownames(Mat)
rownames(d)<-rownames(Mat)
colnames(a)<-c("cox.p","HR","95%CI.low","95%CI.high")
colnames(d)<-c("lnHR","HR","se(lnHR)")
gene.cox[[names(HDP.OVsets)[16]]]<-a
gene.meta[[names(HDP.OVsets)[16]]]<-d

save(gene.cox,file="gene_cox.Rdata")
save(gene.meta,file="gene_meta.Rdata")
save(clinical,file="clinical.Rdata")
save(datasets_nonscale,file="datasets_nonscale.Rdata")
save(datasets_scale,file="datasets_scale.Rdata")
##########get actual numbers of genes related to survival
setwd("G:\\haodapeng-data\\deleteScatter\\think1")
load("gene_cox.Rdata")
threshold<-0.01
actual_genenum<-matrix(sapply(gene.cox,function(x){length(which(x[,1]<threshold))}),nrow=1)
rownames(actual_genenum)<-"genes related to survival"
colnames(actual_genenum)<-names(gene.cox)
write.csv(actual_genenum,"actual_genenum.csv")
diff.coxgene<-lapply(gene.cox,function(x){rownames(x)[which(x[,1]<threshold)]})
save(diff.coxgene,file="diff_coxgene.Rdata")
###########cox gene相似度
source('G:/haodapeng-data/deleteScatter/code/lower2symmetric.R')
source('G:/haodapeng-data/deleteScatter/codes/signatureSimility.R')
load("diff_coxgene.Rdata")
load("datasets_scale.Rdata")
globalExp<-data.frame(GN=rownames(datasets_scale[[1]]),datasets_scale[[1]])
for(i in 2:length(datasets_scale)){
  globalExp<-merge(globalExp,data.frame(GN=rownames(datasets_scale[[i]]),
                                        datasets_scale[[i]]),
                   by.x="GN",by.y="GN",all=T)
}
rownames(globalExp)<-globalExp$GN
globalExp<-globalExp[,-1]
coxGeneSim<-c()
coxGeneSimP<-c()
for(j in 1:15){
  for(i in (j+1):16){
    sa<-diff.coxgene[[i]]
    sb<-diff.coxgene[[j]]
    coxGeneSim<-c(coxGeneSim,signatureSimility(sa,sb,globalExp))
    ramsim<-c()
    for(m in 1:100){
      set.seed(m);ram.a<-sample(rownames(globalExp),length(sa))
      set.seed(m+1);ram.b<-sample(rownames(globalExp),length(sb))
      ramsim[m]<-signatureSimility(ram.a,ram.b,globalExp)
    }
    z_score<-as.numeric(scale(ramsim))
    SD<-sd(ramsim)#
    MEAN<-mean(ramsim)
    real.scale<-as.numeric((coxGeneSim[length(coxGeneSim)]-MEAN)/SD)
    coxGeneSimP<-c(coxGeneSimP,2*pnorm(abs(real.scale),lower.tail=FALSE))
  }
}
coxGeneSim<-lower2symmetric(16,coxGeneSim)
coxGeneSimP<-lower2symmetric(16,coxGeneSimP)
write.csv(coxGeneSim,"coxGeneSim.csv")
write.csv(coxGeneSimP,"coxGeneSimP.csv")
#####得到真实cox基因的变化范围
setwd("G:\\haodapeng-data\\deleteScatter\\think1")
source('G:/haodapeng-data/deleteScatter/code/HGClinicalStatus.R')
source('G:/haodapeng-data/deleteScatter/code/coxAnalysis.R')
library(survival)
load("datasets_scale.Rdata")
load("clinical.Rdata")
gradeStatus<-c(0,rep(1,5),0,rep(1,8))#前15个数据集
result<-HGClinicalStatus(clinical[1:15],gradeStatus)
idx.tumorstage<-result[,1]
idx.tumorgrade<-result[,2]
idx.age<-result[,3]
expNa<-c(rep(0,14),1)

gene.cox.sample<-apply(matrix(1:15,nrow=15),1,function(i){
  phe<-clinical[[i]]
  Mat <- datasets_scale[[i]]
  apply(matrix(1:100,nrow=100),1,function(m){
    set.seed(m);po<-sample(1:nrow(phe),ceiling(nrow(phe)*0.8))
    Mat1<-Mat[,po]
    phe1<-phe[po,]
    result.cox<-coxAnalysis(Mat1,phe1,idx.tumorstage[i],
                            idx.tumorgrade[i],idx.age[i],
                            expNa[i])
    a<-result.cox[[1]]
    length(which(a[,1]<0.01))
  })
})

i=16;
phe<-clinical[[i]]
Mat <- datasets_scale[[i]]
result<-apply(matrix(1:100,nrow=100),1,function(m){
  set.seed(m);po<-sample(1:nrow(phe),ceiling(nrow(phe)*0.8))
  Mat1<-Mat[,po]
  phe1<-phe[po,]
  result1<-apply(matrix(1:nrow(Mat1),nrow=nrow(Mat1)),1,function(j){
    gene.j <- coxph(Surv(donor_survival_time,donor_vital_status == "deceased")
                    ~ Mat1[j,]+as.numeric(tumour_stage)+tumour_grade+
                      donor_age_at_diagnosis, phe1)
    b<-summary(gene.j)
    c(b$coefficients[1,5],b$conf.int[1,c(1,3,4)])
  })
  length(which(result1[1,]<0.01))
})
gene.cox.sample<-cbind(gene.cox.sample,result)


#####cox基因内容重叠度(jaccard)#####
setwd("G:\\haodapeng-data\\deleteScatter\\think1")
load("diff_coxgene.Rdata")
load("datasets_scale.Rdata")
jaccardMatrix<-matrix(nrow=16,ncol=16)
#jaccardMatrix_p<-matrix(nrow=16,ncol=16)
colnames(jaccardMatrix)<-names(diff.coxgene)
rownames(jaccardMatrix)<-names(diff.coxgene)
# colnames(jaccardMatrix_p)<-names(diff.coxgene)
# rownames(jaccardMatrix_p)<-names(diff.coxgene)
for(i in 1:16){
  for(j in 1:16){
    intersect.ij<-length(intersect(diff.coxgene[[i]],diff.coxgene[[j]]))
    union.ij<-length(union(diff.coxgene[[i]],diff.coxgene[[j]]))
    jaccard<-intersect.ij/union.ij
    jaccardMatrix[i,j]<-jaccard
    # randomJaccard.ij<-c()
    # for(m in 1:1000){
    #   set.seed(m);random.i<-sample(rownames(datasets_scale[[i]]),
    #                                length(diff.coxgene[[i]]))
    #   set.seed(m+1);random.j<-sample(rownames(datasets_scale[[j]]),
    #                                length(diff.coxgene[[j]]))
    #   randomInter.ij<-length(intersect(random.i,random.j))
    #   randomUnion.ij<-length(union(random.i,random.j))
    #   randomJaccard.ij<-c(randomJaccard.ij,randomInter.ij/randomUnion.ij)
    # }
    # randomJaccard.ij<-randomJaccard.ij[!is.na(randomJaccard.ij)]
    # z_score<-as.numeric(scale(randomJaccard.ij))
    # SD<-sd(randomJaccard.ij)#
    # MEAN<-mean(randomJaccard.ij)
    # real.scale<-as.numeric((jaccard-MEAN)/SD)
    # jaccardMatrix_p[i,j]<-2*pnorm(abs(real.scale),lower.tail=FALSE)
  }
}
diag(jaccardMatrix)<-NA
#(jaccardMatrix_p)<-NA
write.csv(jaccardMatrix,"jaccardMatrix.csv")
#write.csv(jaccardMatrix_p,"jaccardMatrix_p.csv")
#####随机cox结果#####
setwd("G:\\haodapeng-data\\deleteScatter\\think1")
source('G:/haodapeng-data/deleteScatter/code/HGClinicalStatus.R')
source('G:/haodapeng-data/deleteScatter/code/coxAnalysis.R')
library(GEOquery)
library(survival) 
load("HGSC datasets with OS.rda")
names(HDP.OVsets)[13] <- "Dressman et al.";names(HDP.OVsets)[9]<- "GSE32062"; 
names(HDP.OVsets)[1] <- "E-MTAB-386";names(HDP.OVsets)[16]<- "Patch et al."
gradeStatus<-c(0,rep(1,5),0,rep(1,8))#前15个数据集
#获取各个数据集HG样本的stage,grade,age状态
result<-HGClinicalStatus(HDP.OVsets[1:15],gradeStatus)
idx.tumorstage<-result[,1]
idx.tumorgrade<-result[,2]
idx.age<-result[,3]
expNa<-c(rep(0,14),1)

#获取cox生存相关基因
threshold<-0.01
genenum.random<-list()
for(i in 1:15){
  MyGSE <- HDP.OVsets[[i]];phe <- pData(MyGSE);Mat <- exprs(MyGSE)
  if(gradeStatus[i]==0){
    phe<-phe[which(phe$histological_type=="ser" & phe$sample_type=="tumor" & phe$summarystage=="late"),]
  }else if(gradeStatus[i]==1){
    phe<-phe[which(phe$histological_type=="ser" & phe$sample_type=="tumor" & phe$summarystage=="late" & phe$grade>=2),]
  }
  Mat <- Mat[,rownames(phe)]
  Mat <- t(scale(t(Mat)))
  num.i<-c()
  for (m in 1:100){
    set.seed(m*2+3)
    phe1<-phe[sample(1:nrow(phe),nrow(phe),replace = F),]
    p_value<-coxAnalysis(Mat,phe1,idx.tumorstage[i],
                         idx.tumorgrade[i],idx.age[i],
                         expNa[i])[[1]][,1]
    num.i<-c(num.i,length(which(p_value<threshold)))
  }
  genenum.random[[names(HDP.OVsets)[i]]]<-num.i
}

i=16;MyGSE <- HDP.OVsets[[i]];phe <- pData(MyGSE);Mat <- exprs(MyGSE);
Mat<-Mat[-which(apply(Mat,1,function(x){all(x==x[1])})),]#删掉的是基因
Mat <- t(scale(t(Mat)))
num.i<-c()
for (m in 1:100){
  p_value<-c()
  set.seed(m*2+3)
  phe1<-phe[sample(1:nrow(phe),nrow(phe),replace = F),]
  for (j in 1:nrow(Mat)) {
    gene.j <- coxph(Surv(donor_survival_time,donor_vital_status == "deceased")
                    ~ Mat[j,]+tumour_stage+tumour_grade+
                      donor_age_at_diagnosis, phe1)
    b<-summary(gene.j)$coefficients[1,5]
    p_value <- c(p_value,b)}
  p_value<-p_value[which(!is.na(p_value))]
  num.i<-c(num.i,length(which(p_value<threshold)))
}
genenum.random[[names(HDP.OVsets)[16]]]<-num.i
save(genenum.random,file="genenum_random.Rdata")
num<-sapply(genenum.random,function(x){x})
write.csv(num,"Number of genes in a random situation.csv",row.names = F)
#####实际在随机中的位置######
library(ggplot2)
library(ggpubr)
setwd("G:\\haodapeng-data\\deleteScatter\\think1")
load("genenum_random_m.Rdata")
actual_genenum<-read.csv("actual_genenum.csv",header = T,row.names = 1,as.is=T)
real.random<-c()
for(i in 1:16){
  z_score<-as.numeric(scale(genenum.random[[i]]))
  SD<-sd(genenum.random[[i]])#
  MEAN<-mean(genenum.random[[i]])
  real.scale<-as.numeric((actual_genenum[i]-MEAN)/SD)
  p<-2*pnorm(abs(real.scale),lower.tail=FALSE) 
  real.random<-rbind(real.random,cbind(names(genenum.random)[i],
                                       z_score,real.scale,p))
}
real.random<-data.frame(type=factor(real.random[,1],
                                    levels = names(genenum.random),ordered = T),
                        z_score=as.numeric(real.random[,2]),
                        real=as.numeric(real.random[,3]),
                        p=as.numeric(real.random[,4]),
                        significant=ifelse(as.numeric(real.random[,4])<0.1,
                                           "significant","non-significant"))
bottom<-ggplot(data=real.random,aes(x=type,y=z_score))+
  geom_boxplot(outlier.colour=NA)+
  geom_point(aes(x=type,y=real,colour=significant,fill=significant),
             shape=23,size=2)+
  geom_text(aes(x=type,y=real-.5,label=format(p,digits=2,scientific = T)),size=3.9)+
  labs(x=NULL,y=NULL)+
  scale_color_manual(values=c("black","red"))+
  scale_fill_manual(values=c("black","red"))+
  theme_bw()+
  coord_cartesian(ylim = c(-1.5:3.5))+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))
top<-ggplot(data=real.random,aes(x=type,y=z_score))+
  geom_boxplot(outlier.colour=NA)+
  geom_point(aes(x=type,y=real,colour=significant,fill=significant),
             shape=23,size=2)+
  geom_text(aes(x=type,y=real-.5,label=format(p,digits=2,scientific = T)),size=3.9)+
  labs(x=NULL,y=NULL)+
  scale_color_manual(values=c("black","red"))+
  scale_fill_manual(values=c("black","red"))+
  theme_bw()+
  coord_cartesian(ylim = c(12.5:14))+
  scale_y_continuous(breaks = c(12,13,14))+
  theme(panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
ggarrange(top,bottom,heights=c(1/5,4/5),ncol=1,nrow=2,common.legend = T,legend = "right",align="v")
ggsave("random&actual cox gene.pdf",width = 7,height = 7)

#####cox基因功能富集分析及基因间相似性计算#####
setwd("G:\\haodapeng-data\\deleteScatter\\think1")
library("clusterProfiler")
library("org.Hs.eg.db")
load("diff_coxgene.Rdata")
load("datasets_scale.Rdata")
enrich.result<-lapply(diff.coxgene,function(x){
  h<-as.data.frame(enrichGO(x,
                            OrgDb="org.Hs.eg.db",
                            keytype = "SYMBOL",
                            ont="BP",
                            pvalueCutoff = 1,
                            pAdjustMethod = "fdr",
                            qvalueCutoff = 1,
                            minGSSize = 20,
                            maxGSSize = 5000))
  h1<-h[h[,5]<0.005,]
  return(h1)
})
save(enrich.result,file="enrich_result.Rdata")
#load("enrich_result.Rdata")

###基因间语义相似性###
library("GOSemSim")
allGene<-unique(unlist(sapply(datasets_scale,function(x){rownames(x)})))
d <- godata('org.Hs.eg.db', ont="BP", computeIC=FALSE)
tri_semsim<-c()
tri_semsim_p<-c()
for(i in 14)
{
  for(j  in (1:16)[-14])
  {
    sim<-mgoSim(enrich.result[[i]][,1],enrich.result[[j]][,1],semData=d, measure="Wang")
    tri_semsim<-c(tri_semsim,sim)
    random<-c()
    for(m in 1:100)
    {
      random.ij<-list()
      set.seed(m);random.ij[[1]]<-sample(allGene,length(diff.coxgene[[i]]))
      set.seed(m+1);random.ij[[2]]<-sample(allGene,length(diff.coxgene[[j]]))
      random.enrich<-lapply(random.ij,function(x){
        h<-as.data.frame(enrichGO(x,
                                  OrgDb="org.Hs.eg.db",
                                  keytype = "SYMBOL",
                                  ont="BP",
                                  pvalueCutoff = 1,
                                  pAdjustMethod = "fdr",
                                  qvalueCutoff = 1,
                                  minGSSize = 20,
                                  maxGSSize = 5000))
        h1<-h[h[,5]<0.005,]
        return(h1)
      })
      random.ij<-mgoSim(random.enrich[[1]][,1],random.enrich[[2]][,1],
                        semData=d, measure="Wang")
      random<-c(random,random.ij)
    }
    random<-random[!is.na(random)]
    real.scale<-as.numeric((sim-mean(random))/sd(random))
    p<-2*pnorm(abs(real.scale),lower.tail=FALSE)
    tri_semsim_p<-c(tri_semsim_p,p)
  }
}

semSimilirity<-lower2symmetric(16,tri_semsim)
semSimilirity_p<-lower2symmetric(16,tri_semsim_p)
write.csv(semSimilirity,"semSimilirity.csv")
write.csv(semSimilirity_p,"semSimilirity_p.csv")

#####识别模型基因#####
setwd("G:\\haodapeng-data\\deleteScatter\\think1")
###寻找频率大于2的基因及其meta统计量###
source('G:/haodapeng-data/deleteScatter/code/geneFre.R')
load("datasets_scale.Rdata")
load("gene_meta.Rdata")
trainPo<-which(sapply(datasets_scale,ncol)>100)[-8]
trainSet<-gene.meta[trainPo]
# trainSample<-sum(sapply(datasets_scale[trainPo],ncol))
# testSample<-sum(sapply(datasets_scale[-trainPo],ncol))
gene.metaTrain<-geneFre(trainSet)
save(gene.metaTrain,file="gene_metaTrain.Rdata")
###合并统计量计算###
source('G:/haodapeng-data/deleteScatter/code/meta_summary.R')
library("meta")
library(clusterProfiler)
library(org.Hs.eg.db)
setwd("G:\\haodapeng-data\\deleteScatter\\think1")
load("gene_metaTrain.Rdata")
Train_meta_summary<-meta_summary(gene.metaTrain)
write.csv(Train_meta_summary,"Train_meta_summary.csv")
###风险基因、保护基因获取###
source('G:/haodapeng-data/deleteScatter/code/risk.R')
source('G:/haodapeng-data/deleteScatter/code/protect.R')
Train_meta_summary<-read.csv("Train_meta_summary.csv",
                             header = T,row.names = 1,as.is=T)
colnames(Train_meta_summary)<-NULL
a<-risk(Train_meta_summary)
riskGeneInfor<-a[[1]]
riskGene<-a[[2]]
b<-protect(Train_meta_summary)
protectGeneInfor<-b[[1]]
protectGene<-b[[2]]
write.csv(riskGeneInfor,"riskGeneInfor.csv")
write.table(riskGene,"riskGene.txt",col.names = F,row.names = F)
write.csv(protectGeneInfor,"protectGeneInfor.csv")
write.table(protectGene,"protectGene.txt",col.names = F,row.names = F)
###功能驱动的模型基因识别###
Train_meta_summary<-read.csv("Train_meta_summary.csv",
                             header = T,row.names = 1,as.is=T)
Train_meta_summary1<-Train_meta_summary
colnames(Train_meta_summary1)<-c("Q","df.Q","p.heterogeneity","I2",
                                 "HR","upper","lower","P",
                                 "HR","upper","lower","P")
gene_premodel<-c()
#0为fixed，1为random
for(i in 1:nrow(Train_meta_summary1)){
  if(Train_meta_summary1[i,3]>=0.1)
    gene_premodel<-rbind(gene_premodel,
                         c(unlist(Train_meta_summary1[i,5:8]),0))
  if(Train_meta_summary1[i,3]<0.1)
    gene_premodel<-rbind(gene_premodel,
                         c(unlist(Train_meta_summary1[i,9:12]),1))
}
rownames(gene_premodel)<-rownames(Train_meta_summary1)
gene_premodel<-cbind(gene_premodel,
                          p.adjust(gene_premodel[,4],method = "fdr"))
colnames(gene_premodel)<-c("HR","upper","lower","P","model","p.adjust")
gene_enrich<-gene_premodel[which(gene_premodel[,6]<0.01),]
write.csv(gene_enrich,"gene_enrich.csv")
geneNamePreModel<-rownames(gene_enrich)
write.table(geneNamePreModel,"geneNamePreModel.txt",quote = F,col.names = F,row.names = F)
#geneNamePreModel<-as.character(unlist(read.table("geneNamePreModel.txt",as.is=T)))
h<-enrichGO(geneNamePreModel,
                      OrgDb="org.Hs.eg.db",
                      keytype = "SYMBOL",
                      ont="BP",
                      pvalueCutoff = 1,
                      pAdjustMethod = "fdr",
                      minGSSize = 500,
                      maxGSSize = 2000)
enrich_meta<-h[h[,7]<0.001,]
write.csv(enrich_meta,"enrich_meta.csv")
enrich_gene<-unique(unlist(strsplit(enrich_meta$geneID,"/")))
write.table(enrich_gene,"enrich_gene.txt",quote = F,row.names = F,col.names = F)
#enrich_meta<-read.csv("enrich_meta.csv",header=T,row.names=1,as.is=T)
enrich_meta$Significant1<-(-log(enrich_meta$qvalue,base=10))
enrich_meta<-as.data.frame(enrich_meta[1:20,])
ggplot(data = enrich_meta)+
  geom_bar(aes(x=reorder(Description,Significant1),y=enrich_meta[["Significant1"]],
               fill=enrich_meta[["Significant1"]]),stat = 'identity',width=0.9)+
  geom_hline(yintercept = -log(0.05,base=10),linetype=2,size=0.6)+
  coord_flip()+
  scale_fill_gradient(low = "#92B4D7",high = "#2171B5")+
  xlab("")+
  ylab(expression(paste(-log[10],"(q value)")))+
  theme_bw()+
  theme(
    axis.text.x = element_text(color="black",size = 12),
    axis.text.y = element_text(color="black",size = 12),
    axis.title.x = element_text(color="black",size = 16),
    legend.text = element_text(color="black",size = 12),
    legend.title =element_text(color="black",size = 12),
    panel.border = element_blank(),panel.grid=element_blank(),
    axis.line = element_line(colour = "black")
  )+
  labs(fill="significant\n")
ggsave("功能富集图Q.pdf",width=7,height=5)

enrich_meta$Significant2<-(-log(enrich_meta$pvalue))
ggplot(data = enrich_meta)+
  geom_bar(aes(x=reorder(Description,Significant2),
               y=enrich_meta[["Significant2"]],fill=enrich_meta[["Significant2"]]),stat = 'identity')+
  coord_flip()+
  scale_fill_gradient(low = "#FFBD00",high = "#FF0000")+
  xlab("")+
  ylab("-logPvalue")+
  theme_bw()+
  theme(
    axis.text.x = element_text(color="black",size = rel(1.5)),
    axis.text.y = element_text(color="black",size = rel(1.3)),
    axis.title.x = element_text(color="black",size = rel(1.6)),
    legend.text = element_text(color="black",size = rel(1.0)),
    legend.title =element_text(color="black",size = rel(1.1))
  )+
  labs(fill="significant\n")
ggsave("功能富集图P.pdf",width=10,height=6)


###把模型基因划为保护，风险两类###
enrich_gene<-as.character(unlist(read.table("enrich_gene.txt",as.is=T)))
protectGene<-as.character(unlist(read.table("protectGene.txt",as.is=T)))
riskGene<-as.character(unlist(read.table("riskGene.txt",as.is=T)))
protectModelGene<-c()
riskModelGene<-c()
for(i in 1:length(enrich_gene)){
  if(enrich_gene[i]%in%protectGene){
    protectModelGene<-c(protectModelGene,enrich_gene[i])
  }else if(enrich_gene[i]%in%riskGene){
    riskModelGene<-c(riskModelGene,enrich_gene[i])
  }
}
write.table(protectModelGene,"protectModelGene.txt",quote = F,row.names = F,col.names = F)
write.table(riskModelGene,"riskModelGene.txt",quote = F,row.names = F,col.names = F)
###为模型基因分配权重###
protectModelGene<-as.character(unlist(read.table("protectModelGene.txt",as.is=T)))
riskModelGene<-as.character(unlist(read.table("riskModelGene.txt",as.is=T)))
protectGeneInfor<-read.csv("protectGeneInfor.csv",as.is=T,header = T,row.names = 1)
riskGeneInfor<-read.csv("riskGeneInfor.csv",as.is=T,header = T,row.names = 1)
protectWeight<-abs(log(protectGeneInfor[protectModelGene,1]))
riskWeight<-abs(log(riskGeneInfor[riskModelGene,1]))
#protectWeight<-(-log(protectGeneInfor[protectModelGene,4]))
write.table(protectWeight,"protectWeight.txt",quote = F,row.names = F,col.names = F)
write.table(riskWeight,"riskWeight.txt",quote = F,row.names = F,col.names = F)


####统计####
gene_enrich<-read.csv("gene_enrich.csv",header=T,row.names = 1,as.is=T)
protectModelGene<-as.character(unlist(read.table("protectModelGene.txt",as.is=T)))
riskModelGene<-as.character(unlist(read.table("riskModelGene.txt",as.is=T)))
gene_enrich_68<-gene_enrich[c(protectModelGene,riskModelGene),]
gene_enrich_68$class<-c(rep("Protect",43),rep("Risk",25))
gene_enrich_68$Weight<-abs(log(gene_enrich_68$HR))
a<-gene_enrich_68[,-5]
b<-data.frame(HR=format(a$HR,digits=3),
              CI=paste(format(a$lower,digits=3),
                       format(a$upper,digits=3),sep="-"),
              Pvalue=format(a$P,digits=2,scientific = 3),
              Padjust=format(a$p.adjust,digits=2),
              Class=a$class,
              Weight=format(a$Weight,digits=3))
rownames(b)<-rownames(gene_enrich_68)
write.csv(b,"gene_enrich_68.csv")
#####模型应用#####
source('G:/haodapeng-data/deleteScatter/code/double-weighted_ssGSEA.R')
source('G:/haodapeng-data/deleteScatter/code/ssgseaScore.R')
setwd("G:\\haodapeng-data\\deleteScatter\\think1")
load("datasets_nonscale.Rdata")
load("datasets_scale.Rdata")
load("TCGA_seq.Rdata")
protectModelGene<-as.character(unlist(read.table("protectModelGene.txt",as.is=T)))
riskModelGene<-as.character(unlist(read.table("riskModelGene.txt",as.is=T)))
protectWeight<-as.numeric(unlist(read.table("protectWeight.txt",as.is=T)))
riskWeight<-as.numeric(unlist(read.table("riskWeight.txt",as.is=T)))
trainPo<-which(sapply(datasets_scale,ncol)>100)[-8]
testPo<-which(sapply(datasets_scale,ncol)<100)[-8]
independPo<-c(14,16)
ssgsea.out<-ssgseaScore(datasets_nonscale[trainPo],w1=1.5,gProtect = protectModelGene,
                        gRisk = riskModelGene,wProtect = protectWeight,
                        wRisk = riskWeight)
ssgsea.out.independ<-ssgseaScore(datasets_nonscale[independPo],w1=1.5,gProtect = protectModelGene,
                                 gRisk = riskModelGene,wProtect = protectWeight,
                                 wRisk = riskWeight)
ssgsea.out.test<-ssgseaScore(datasets_nonscale[testPo],w1=1.5,gProtect = protectModelGene,
                             gRisk = riskModelGene,wProtect = protectWeight,
                             wRisk = riskWeight)
ssgsea.out.all<-ssgseaScore(datasets_nonscale,w1=1.5,gProtect = protectModelGene,
                            gRisk = riskModelGene,wProtect = protectWeight,
                            wRisk = riskWeight)
# ssgsea.out.tcga<-ssgseaScore(list(TCGA.seq[[1]]),w1=1.5,gProtect = protectModelGene,
#                              gRisk = riskModelGene,wProtect = protectWeight,
#                              wRisk = riskWeight)
# ssgsea.out.patch<-ssgseaScore(datasets_nonscale[-trainPo][9],w1=1.5,gProtect = protectModelGene,
#                               gRisk = riskModelGene,wProtect = protectWeight,
#                               wRisk = riskWeight)
# ssgsea.out.tcga.cel<-ssgseaScore(datasets_nonscale[14],w1=1.5,gProtect = protectModelGene,
#                                  gRisk = riskModelGene,wProtect = protectWeight,
#                                  wRisk = riskWeight)
save(ssgsea.out,file="ssgsea_out.Rdata")
save(ssgsea.out.test,file="ssgsea_out_test.Rdata")
save(ssgsea.out.independ,file="ssgsea.out.independ.Rdata")
save(ssgsea.out.all,file="ssgsea_out_all.Rdata")
# save(ssgsea.out.tcga,file="ssgsea.out.tcga.Rdata")
# save(ssgsea.out.patch,file="ssgsea.out.patch.Rdata")
# save(ssgsea.out.tcga.cel,file="ssgsea.out.tcga.cel.Rdata")
#####模型评估#####
source('G:/haodapeng-data/deleteScatter/code/signatureScore.R')
source('G:/haodapeng-data/deleteScatter/code/getTotalGene.R')
source('G:/haodapeng-data/deleteScatter/code/ROC_CI_global.R')
source('G:/haodapeng-data/deleteScatter/code/globalEval.R')
source('G:/haodapeng-data/deleteScatter/code/HGClinicalStatus.R')
source('G:/haodapeng-data/deleteScatter/code/plotSurvivalRoc.R')
source('G:/haodapeng-data/deleteScatter/code/sstselfEvalUMContinue.R')
source('G:/haodapeng-data/deleteScatter/code/selfEval.R')
source('G:/haodapeng-data/deleteScatter/code/sstselfEvalUM.R')
source('G:/haodapeng-data/deleteScatter/code/unicoxForestData.R')
source('G:/haodapeng-data/deleteScatter/code/MulticoxForestData.R')
source('G:/haodapeng-data/deleteScatter/code/sstselfEval.R')
source('G:/haodapeng-data/deleteScatter/code/groupQ.R')
source('G:/haodapeng-data/deleteScatter/code/createTable.R')
source('G:/haodapeng-data/deleteScatter/code/ROC_CI.R')
source('G:/haodapeng-data/deleteScatter/think1/combineCoxResult.R')
source('G:/haodapeng-data/deleteScatter/code/getSurvivalMatrix.R')
source('G:/haodapeng-data/deleteScatter/code/plotSurvivalCurve.R')
source('G:/haodapeng-data/deleteScatter/code/getSurvivalMatrix1.R')
source('G:/haodapeng-data/deleteScatter/code/group.R')
source('G:/haodapeng-data/deleteScatter/code/Digits4.R')
source('G:/haodapeng-data/deleteScatter/code/getGlobalExp.R')
source('G:/haodapeng-data/deleteScatter/code/getGlobalScore.R')
source('G:/haodapeng-data/deleteScatter/code/getHeatmapExp.R')
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
setwd("G:\\haodapeng-data\\deleteScatter\\think1")
load("model_method.RData")
load("model_coefs1.Rdata")
load("datasets_nonscale.Rdata")
load("datasets_scale.Rdata")
load("ssgsea_out.Rdata")
load("ssgsea_out_test.Rdata")
load("ssgsea_out_all.Rdata")
load("ssgsea.out.independ.Rdata")
# load("ssgsea.out.tcga.Rdata")
# load("ssgsea.out.tcga.cel.Rdata")
# load("ssgsea.out.patch.Rdata")
load("clinical.Rdata")
protectModelGene<-as.character(unlist(read.table("protectModelGene.txt",as.is=T)))
riskModelGene<-as.character(unlist(read.table("riskModelGene.txt",as.is=T)))
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
# score.test.tcga<-structure(as.numeric(unlist(sapply(ssgsea.out.tcga,function(x){x[3,]}))),
#                            .Names=as.character(unlist(sapply(ssgsea.out.tcga,colnames))))
score.test<-structure(as.numeric(unlist(sapply(ssgsea.out.test,function(x){x[3,]}))),
                      .Names=as.character(unlist(sapply(ssgsea.out.test,colnames))))
score.all<-structure(as.numeric(unlist(sapply(ssgsea.out.all,function(x){x[3,]}))),
                      .Names=as.character(unlist(sapply(ssgsea.out.all,colnames))))
score.test.independ<-structure(as.numeric(unlist(sapply(ssgsea.out.independ,function(x){x[3,]}))),
                               .Names=as.character(unlist(sapply(ssgsea.out.independ,colnames))))
# score.test.patch<-structure(as.numeric(unlist(sapply(ssgsea.out.patch,function(x){x[3,]}))),
#                             .Names=as.character(unlist(sapply(ssgsea.out.patch,colnames))))
# score.test.tcga.cel<-structure(as.numeric(unlist(sapply(ssgsea.out.tcga.cel,function(x){x[3,]}))),
#                            .Names=as.character(unlist(sapply(ssgsea.out.tcga.cel,colnames))))
#####1.训练集生存数据####
clinicalTrain<-clinical[trainPo]
clinicalMatrix<-getSurvivalMatrix(clinicalTrain,"global","Patch et al.")
survivalMatrix<-as.data.frame(cbind(clinicalMatrix,score,1))
write.csv(survivalMatrix,"survivalMatrix.csv")

#####2.分数分布表#####
a<-sapply(ssgsea.out,function(x){
  np<-length(which(x[3,]>medianScore))
  nn<-length(which(x[3,]<medianScore))
  return(c(np,nn))
})
rownames(a)<-c("score>median score","score<median score")

#####3.分数分布曲线图#####
######训练集#####
survivalMatrix<-read.csv("survivalMatrix.csv",header = T,row.names = 1,as.is = T)
survivalMatrix_sort<-survivalMatrix[order(survivalMatrix[,3]),]
survivalMatrix_sort<-cbind(survivalMatrix_sort,
                           seq(from=1,to=nrow(survivalMatrix_sort),by=1))
colnames(survivalMatrix_sort)[c(2,5)]<-c("SurvivalStatus","order")
survivalMatrix_sort$score<-survivalMatrix_sort$score-medianScore

poL0<-which(survivalMatrix_sort$score>0)
poS0<-which(survivalMatrix_sort$score<0)
ma<-matrix(c(length(which((survivalMatrix_sort$score>0)&(survivalMatrix_sort$SurvivalStatus=="Alive"))),
             length(which((survivalMatrix_sort$score<0)&(survivalMatrix_sort$SurvivalStatus=="Alive"))),
             length(which((survivalMatrix_sort$score>0)&(survivalMatrix_sort$SurvivalStatus=="Dead"))),
             length(which((survivalMatrix_sort$score<0)&(survivalMatrix_sort$SurvivalStatus=="Dead")))),nrow = 2)
rownames(ma)<-c("High Risk","Low Risk")
colnames(ma)<-c("Alive","Dead")
chisq.test(ma)

ggplot(survivalMatrix_sort,aes(order,score))+
  geom_bar(stat = "identity",aes(colour=SurvivalStatus,fill=SurvivalStatus))+
  scale_color_manual(values=c("#F5BB2C","#033994"))+
  scale_fill_manual(values =c("#F5BB2C","#033994"))+
  theme_classic()+
  scale_y_continuous(breaks  = c(-12360,-2688.0,11330)-medianScore,
                     labels = c("0","50th","100th"))+
  labs(x="",y="Percentile of FIPI",title="Training")+
  geom_vline(xintercept = 595,linetype="dotted")+
  theme(aspect.ratio =7/10,legend.position = c(0.01,1),
        legend.justification = c(0.01,1),
        plot.title = element_text(hjust = 0.5,size=15),
        legend.text = element_text(size=12),
        legend.title = element_text(size=14),
        axis.text.x = element_text(color = "black",size=14),
        axis.text.y=element_text(color="black",size=14),
        axis.title.y=element_text(color="black",size=16))+
  annotate("text",length(score)/2-400,-6000,label="Low-Risk",color="#01939A",size=5)+
  annotate("text",length(score)/2+400,6100,label="High-Risk",color="#E73018",size=5)
ggsave("Training Score distribution.pdf",width=6,height = 4)

ggplot(survivalMatrix_sort,aes(order,time))+
  geom_point(aes(color=SurvivalStatus))+
  labs(x="",y="Survival Time (Day)")+
  geom_vline(xintercept = 1686/2)+
  theme_classic()+
  theme(aspect.ratio = 9/10,
        legend.position = c(0.01,1),
        legend.justification = c(0.01,1),
        legend.text = element_text(size=12),
        legend.title = element_text(size=14),
        axis.text.x = element_text(color = "black",size=14),
        axis.text.y=element_text(color="black",size=14),
        axis.title.y=element_text(color="black",size=16))+
  scale_colour_manual(values =c("#FFAC40",rgb(0,0,0)))
ggsave("分布散点图.pdf",width=6,height=6)

Data<-getHeatmapExp(gene,datasets_nonscale[trainPo])
Data1<-Data[,rownames(survivalMatrix_sort)]
annotation_row<-data.frame(Gene=c(rep("Protect",25),rep("Risk",43)))
rownames(annotation_row)<-gene
annotation_col<-data.frame(Sample=c(rep("High Risk",595),
                                   rep("Low Risk",595)))
rownames(annotation_col)<-colnames(Data1)
pheatmap(Data1,cluster_rows = F,cluster_cols = F,scale="row",
         show_colnames = F,annotation_row = annotation_row,
         annotation_col = annotation_col,na_col="green")
######测试集#####
clinicalTest<-clinical[testPo]
clinicalMatrix<-getSurvivalMatrix(clinicalTest,"global","Patch et al.")
# group.test<-c()
# for(i in 1:length(score.test)){
#   group.test.i<-group(score.test[i],type="0",medianScore)
#   group.test<-c(group.test,group.test.i)
# }
# survivalMatrix<-data.frame(time=as.numeric(clinicalMatrix[,1]),
#                            status=as.character(clinicalMatrix[,2]),
#                            score=-as.numeric(score.test),
#                            group=group.test)
survivalMatrix<-data.frame(time=as.numeric(clinicalMatrix[,1]),
                           status=as.character(clinicalMatrix[,2]),
                           score=as.numeric(score.test))
#names(survivalMatrix)<-c("time","status","score","group")
#survivalMatrix$status<-ifelse(survivalMatrix$status=="Alive",0,1)
#survivalMatrix$group<-ifelse(survivalMatrix$group=="negative","High Risk","Low Risk")
survivalMatrix_sort<-survivalMatrix[order(survivalMatrix[,3]),]
survivalMatrix_sort<-cbind(survivalMatrix_sort,
                           seq(from=1,to=nrow(survivalMatrix_sort),by=1))
colnames(survivalMatrix_sort)[c(2,4)]<-c("SurvivalStatus","order")
# survivalMatrix_sort$SurvivalStatus<-factor(survivalMatrix_sort$SurvivalStatus)
survivalMatrix_sort$score<-survivalMatrix_sort$score-medianScore#medianScore

poL0<-which(survivalMatrix_sort$score>0)
poS0<-which(survivalMatrix_sort$score<0)
ma<-matrix(c(length(which((survivalMatrix_sort$score>0)&(survivalMatrix_sort$SurvivalStatus=="Alive"))),
             length(which((survivalMatrix_sort$score<0)&(survivalMatrix_sort$SurvivalStatus=="Alive"))),
             length(which((survivalMatrix_sort$score>0)&(survivalMatrix_sort$SurvivalStatus=="Dead"))),
             length(which((survivalMatrix_sort$score<0)&(survivalMatrix_sort$SurvivalStatus=="Dead")))),nrow = 2)
rownames(ma)<-c("High Risk","Low Risk")
colnames(ma)<-c("Alive","Dead")
chisq.test(ma)


ggplot(survivalMatrix_sort,aes(order,score))+
  geom_bar(stat = "identity",aes(colour=SurvivalStatus,fill=SurvivalStatus))+
  scale_color_manual(values=c("#F5BB2C","#033994"))+
  scale_fill_manual(values =c("#F5BB2C","#033994"))+
  theme_classic()+
  labs(x="",y="",title="Test")+
  geom_vline(xintercept = 195,linetype="dotted")+
  theme(aspect.ratio =7/10,legend.position = c(0.01,1),
        legend.justification = c(0.01,1),
        plot.title = element_text(hjust = 0.5,size=15),
        legend.text = element_text(size=12),
        legend.title = element_text(size=14),
        axis.text.x = element_text(color = "black",size=14),
        #axis.text.y=element_blank(),
        axis.title.y=element_text(color="black",size=16))+
  annotate("text",100,-6000,label="Low-Risk",color="#01939A",size=5)+
  annotate("text",350,6100,label="High-Risk",color="#E73018",size=5)
ggsave("Test Score distribution.pdf",width=6,height = 4)
######独立测试集######
clinicalIndepend<-clinical[independPo]#TCGA、patch
clinicalMatrix11<-getSurvivalMatrix1(clinicalIndepend[1],"global","Patch et al.")
clinicalMatrix22<-getSurvivalMatrix(clinicalIndepend[2],"global","Patch et al.")
clinicalMatrix<-rbind(clinicalMatrix11,clinicalMatrix22)
group.test.independ<-c()
for(i in 1:length(score.test.independ)){
  group.test.i<-group(score.test.independ[i],type="0",medianScore)
  group.test.independ<-c(group.test.independ,group.test.i)
}
survivalMatrix<-data.frame(time=as.numeric(clinicalMatrix[,1]),
                           status=as.character(clinicalMatrix[,2]),
                           score=as.numeric(score.test.independ),
                           group=group.test.independ)
names(survivalMatrix)<-c("time","status","score","group")
#survivalMatrix$status<-ifelse(survivalMatrix$status=="Alive",0,1)
survivalMatrix$group<-ifelse(survivalMatrix$group=="negative","High Risk","Low Risk")
survivalMatrix_sort<-survivalMatrix[order(survivalMatrix[,3]),]
survivalMatrix_sort<-cbind(survivalMatrix_sort,
                           seq(from=1,to=nrow(survivalMatrix_sort),by=1))
colnames(survivalMatrix_sort)[c(2,5)]<-c("SurvivalStatus","order")
survivalMatrix_sort$SurvivalStatus<-factor(survivalMatrix_sort$SurvivalStatus)
survivalMatrix_sort$score<-survivalMatrix_sort$score-medianScore
ggplot(survivalMatrix_sort,aes(order,score))+
  geom_bar(stat = "identity",aes(colour=SurvivalStatus,fill=SurvivalStatus))+
  scale_color_manual(values=c("#F5BB2C","#033994"))+
  scale_fill_manual(values =c("#F5BB2C","#033994"))+
  theme_classic()+
  labs(x="",y="",title="Independent Test")+
  geom_vline(xintercept = 159,linetype="dotted")+
  theme(aspect.ratio =7/10,legend.position = c(0.01,1),
        legend.justification = c(0.01,1),
        legend.text = element_text(size=12),
        plot.title = element_text(hjust = 0.5,size=15),
        legend.title = element_text(size=14),
        axis.text.x = element_text(color = "black",size=14),
        axis.text.y=element_blank(),
        axis.title.y=element_text(color="black",size=16))+
  annotate("text",100,-4000,label="Low-Risk",color="#01939A",size=5)+
  annotate("text",250,4000,label="High-Risk",color="#E73018",size=5)
ggsave("Independent Test Score distribution.pdf",width=6,height = 4)
#####4.生存曲线#####
######训练集整体######
survivalMatrix<-read.csv("survivalMatrix.csv",header = T,row.names = 1,as.is = T)
for(i in 1:nrow(survivalMatrix)){
  survivalMatrix[i,4]<-group(survivalMatrix[i,3],type="0",medianScore)
}
names(survivalMatrix)<-c("time","status","score","group")
survivalMatrix$status<-ifelse(survivalMatrix$status=="Alive",0,1)
survivalMatrix$group<-ifelse(survivalMatrix$group=="negative","High Risk","Low Risk")
eee<-table(survivalMatrix$status,survivalMatrix$group)
chisq.test(eee)

fit<-survfit(Surv(time,status)~group,data=survivalMatrix)
diff<-survdiff(Surv(time,status)~group,data=survivalMatrix)
pdf("Train Global survival curve.pdf",height=6,width=6)
ggsurvplot(fit, data = survivalMatrix, title="Training dataset(n = 1190)",
           risk.table = TRUE,risk.table.col = "strata",pval = 1e-10,
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
dev.off()

#四分
survivalMatrix<-read.csv("survivalMatrix.csv",header = T,row.names = 1,as.is = T)
for(i in 1:nrow(survivalMatrix)){
  survivalMatrix[i,4]<-groupQ(survivalMatrix[i,3],Q1score,medianScore,Q3score)
}
names(survivalMatrix)<-c("time","status","score","group")
survivalMatrix$status<-ifelse(survivalMatrix$status=="Alive",0,1)
survivalMatrix$group<-factor(survivalMatrix$group,levels = c("Class one","Class two",
                                                             "Class three","Class four"))
eee<-table(survivalMatrix$status,survivalMatrix$group)
chisq.test(eee)
fit<-survfit(Surv(time,status)~group,data=survivalMatrix)
diff<-survdiff(Surv(time,status)~group,data=survivalMatrix)#2.22e-16
pdf("Train Global survival curve.pdf",height=6,width=6)
ggsurvplot(fit, data = survivalMatrix, title="Training dataset(n = 1190)",
           risk.table = TRUE,risk.table.col = "strata",pval = T,
           risk.table.height = 0.2,xlab="Time(Days)",legend.title="",
           #legend.labs=c("High Risk","Low Risk"),
           risk.table.title="No. at risk",
           #palette = c("#E7201C","#0A499E"),
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
dev.off()
# 训练集单独
clinicalTrain<-clinical[trainPo]
sink("information of Train seperate survival curve.txt",append=T)
pdf("Train seperate survival curve.pdf",height=6,width=6)
for(i in 1:length(ssgsea.out)){
  survivalMatrix.i<-getSurvivalMatrix(clinicalTrain[[i]],
                                      "seperate","Patch et al.",
                                      names(clinicalTrain)[i])
  survivalMatrix.ii<-cbind(survivalMatrix.i,ssgsea.out[[i]][3,],1)
  for(j in 1:nrow(survivalMatrix.ii)){
    survivalMatrix.ii[j,4]<-group(as.numeric(survivalMatrix.ii[j,3]),
                                  type="0",medianScore)
  }
  survivalMatrix.iii<-data.frame(time=as.numeric(survivalMatrix.ii[,1]),
                                 status=as.character(survivalMatrix.ii[,2]),
                                 score=as.numeric(survivalMatrix.ii[,3]),
                                 group=as.character(survivalMatrix.ii[,4]))
  survivalMatrix.iii$status<-ifelse(survivalMatrix.iii$status=="Alive",0,1)
  survivalMatrix.iii$group<-ifelse(survivalMatrix.iii$group=="negative",
                                   "Low Risk","High Risk")
  fit<-survfit(Surv(time,as.numeric(status))~group,
               data=survivalMatrix.iii)
  diff<-survdiff(Surv(time,as.numeric(status))~group,
               data=survivalMatrix.iii)
  print(names(ssgsea.out)[i])
  print(fit)
  print(diff)
  ggsurvplot(fit, data = survivalMatrix.iii, title=paste(names(ssgsea.out)[i],"(n = ",
                                                         sapply(ssgsea.out,ncol)[i],")",sep=""),
             risk.table = TRUE,risk.table.col = "strata",pval = T,
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
}
dev.off()
sink()
###TCGA-cel
clinicalMatrix<-getSurvivalMatrix(clinical[14],"global","Patch et al.")
colnames(clinicalMatrix)<-c("time","status")
for(i in 1:nrow(clinicalMatrix)){
  if(clinicalMatrix[i,2]=="Dead")
    clinicalMatrix[i,2]<-1
  if(clinicalMatrix[i,2]=="Alive")
    clinicalMatrix[i,2]<-"0"
}
group.test<-c()
for(i in 1:length(score.test.tcga.cel)){
  group.test.i<-group(score.test.tcga.cel[i],type="0",medianScore)
  group.test<-c(group.test,group.test.i)
}
survivalMatrix<-data.frame(time=as.numeric(clinicalMatrix[,1]),
                           status=as.character(clinicalMatrix[,2]),
                           score=score.test.tcga.cel,
                           group=group.test)
names(survivalMatrix)<-c("time","status","score","group")
survivalMatrix$group<-ifelse(survivalMatrix$group=="negative","High Risk","Low Risk")
summary(coxph(Surv(time,as.numeric(status))~ score,survivalMatrix))$coefficients[-4]
fit<-survfit(Surv(time,as.numeric(status))~group,data=survivalMatrix)
diff<-survdiff(Surv(time,as.numeric(status))~group,data=survivalMatrix)
pdf("Test TCGAcel survival curve.pdf",height=6,width=6)
ggsurvplot(fit, data = survivalMatrix, title="Test dataset(n = 495)",
           risk.table = TRUE,risk.table.col = "strata",
           pval =  0.0343,
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
dev.off() 
######测试集整体######
clinicalTest<-clinical[testPo]
clinicalMatrix<-getSurvivalMatrix(clinicalTest,"global","Patch et al.")
group.test<-c()
for(i in 1:length(score.test)){
  group.test.i<-group(score.test[i],type="0",medianScore)
  group.test<-c(group.test,group.test.i)
}
survivalMatrix<-data.frame(time=as.numeric(clinicalMatrix[,1]),
                           status=as.character(clinicalMatrix[,2]),
                           score=as.numeric(score.test),
                           group=group.test)
names(survivalMatrix)<-c("time","status","score","group")
survivalMatrix$status<-ifelse(survivalMatrix$status=="Alive",0,1)
survivalMatrix$group<-ifelse(survivalMatrix$group=="negative","High Risk","Low Risk")
eee<-table(survivalMatrix$status,survivalMatrix$group)
chisq.test(eee)
fit<-survfit(Surv(time,as.numeric(status))~group,data=survivalMatrix)
diff<-survdiff(Surv(time,as.numeric(status))~group,data=survivalMatrix)
pdf("Test Global survival curve.pdf",height=6,width=6)
ggsurvplot(fit, data = survivalMatrix, title="Test dataset(n = 438)",
           risk.table = TRUE,risk.table.col = "strata",
           pval =  0.000364,
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
dev.off()

#四分
clinicalTest<-clinical[testPo]
clinicalMatrix<-getSurvivalMatrix(clinicalTest,"global","Patch et al.")
group.test<-c()
for(i in 1:length(score.test)){
  group.test.i<-groupQ(score.test[i],Q1score,medianScore,Q3score)
  group.test<-c(group.test,group.test.i)
}
survivalMatrix<-data.frame(time=as.numeric(clinicalMatrix[,1]),
                           status=as.character(clinicalMatrix[,2]),
                           score=as.numeric(score.test),
                           group=group.test)
names(survivalMatrix)<-c("time","status","score","group")
survivalMatrix$status<-ifelse(survivalMatrix$status=="Alive",0,1)
survivalMatrix$group<-factor(survivalMatrix$group,levels = c("Class one","Class two",
                                                             "Class three","Class four"))
eee<-table(survivalMatrix$status,survivalMatrix$group)
chisq.test(eee)
fit<-survfit(Surv(time,as.numeric(status))~group,data=survivalMatrix)
diff<-survdiff(Surv(time,as.numeric(status))~group,data=survivalMatrix)
pdf("Test Global survival curve.pdf",height=6,width=6)
ggsurvplot(fit, data = survivalMatrix, title="Test dataset(n = 438)",
           risk.table = TRUE,risk.table.col = "strata",
           pval =  T,
           risk.table.height = 0.2,xlab="Time(Days)",legend.title="",
           #legend.labs=c("High Risk","Low Risk"),
           risk.table.title="No. at risk",
           #palette = c("#E7201C","#0A499E"),
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
dev.off()
#######independ######
clinicalIndepend<-clinical[independPo]#TCGA、patch
clinicalMatrix11<-getSurvivalMatrix1(clinicalIndepend[1],"global","Patch et al.")
clinicalMatrix22<-getSurvivalMatrix(clinicalIndepend[2],"global","Patch et al.")
clinicalMatrix<-rbind(clinicalMatrix11,clinicalMatrix22)
group.test.independ<-c()
for(i in 1:length(score.test.independ)){
  group.test.i<-group(score.test.independ[i],type="0",medianScore)
  group.test.independ<-c(group.test.independ,group.test.i)
}
survivalMatrix<-data.frame(time=as.numeric(clinicalMatrix[,1]),
                           status=as.character(clinicalMatrix[,2]),
                           score=as.numeric(score.test.independ),
                           group=group.test.independ)
names(survivalMatrix)<-c("time","status","score","group")
survivalMatrix$group<-ifelse(survivalMatrix$group=="negative","High Risk","Low Risk")
eee<-table(survivalMatrix$status,survivalMatrix$group)
chisq.test(eee)
fit<-survfit(Surv(time,as.numeric(status))~group,data=survivalMatrix)
diff<-survdiff(Surv(time,as.numeric(status))~group,data=survivalMatrix)
pdf("independent Test survival curve.pdf",height=6,width=6)
ggsurvplot(fit, data = survivalMatrix, title=" Independent Test Dataset(n = 354)",
           risk.table = TRUE,risk.table.col = "strata",
           pval =  0.00359,
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
dev.off()

#四分
clinicalIndepend<-clinical[independPo]#TCGA、patch
clinicalMatrix11<-getSurvivalMatrix1(clinicalIndepend[1],"global","Patch et al.")
clinicalMatrix22<-getSurvivalMatrix(clinicalIndepend[2],"global","Patch et al.")
clinicalMatrix<-rbind(clinicalMatrix11,clinicalMatrix22)
group.test.independ<-c()
for(i in 1:length(score.test.independ)){
  group.test.i<-groupQ(score.test.independ[i],Q1score,medianScore,Q3score)
  group.test.independ<-c(group.test.independ,group.test.i)
}
survivalMatrix<-data.frame(time=as.numeric(clinicalMatrix[,1]),
                           status=as.character(clinicalMatrix[,2]),
                           score=as.numeric(score.test.independ),
                           group=group.test.independ)
names(survivalMatrix)<-c("time","status","score","group")
survivalMatrix$group<-factor(survivalMatrix$group,levels = c("Class one","Class two",
                                                             "Class three","Class four"))
eee<-table(survivalMatrix$status,survivalMatrix$group)
chisq.test(eee)
fit<-survfit(Surv(time,as.numeric(status))~group,data=survivalMatrix)
diff<-survdiff(Surv(time,as.numeric(status))~group,data=survivalMatrix)
pdf("independent Test survival curve.pdf",height=6,width=6)
ggsurvplot(fit, data = survivalMatrix, title=" Independent Test Dataset(n = 354)",
           risk.table = TRUE,risk.table.col = "strata",
           pval = T,
           risk.table.height = 0.2,xlab="Time(Days)",legend.title="",
           #legend.labs=c("High Risk","Low Risk"),
           risk.table.title="No. at risk",
           #palette = c("#E7201C","#0A499E"),
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
dev.off()
###测试集单独
clinicalTest<-clinical[-trainPo]
pdf("Test seperate survival curve.pdf",height=7,width=7)
for(i in 1:length(ssgsea.out.test)){
  survivalMatrix.i<-getSurvivalMatrix(clinicalTest[[i]],
                                      "seperate","Patch et al.",
                                      names(clinicalTest)[i])
  survivalMatrix.ii<-cbind(survivalMatrix.i,ssgsea.out.test[[i]][4,],1)
  for(j in 1:nrow(survivalMatrix.ii)){
    survivalMatrix.ii[j,4]<-group(as.numeric(survivalMatrix.ii[j,3]),
                                  type="0",summary(score)[-4])
  }
  survivalMatrix.iii<-data.frame(time=as.numeric(survivalMatrix.ii[,1]),
                                 status=as.character(survivalMatrix.ii[,2]),
                                 score=as.numeric(survivalMatrix.ii[,3]),
                                 group=as.character(survivalMatrix.ii[,4]))
  fit<-survfit(Surv(time,as.numeric(status))~group,
               data=survivalMatrix.iii)
  ggsurvplot(fit, data = survivalMatrix.iii, size = 1,
             risk.table = TRUE,risk.table.col = "strata",pval = TRUE,
             risk.table.height = 0.2,
             ggtheme = theme_bw())
}
dev.off()
#####5.训整+测1+测2自我评估#####
#####按临床分亚组####
#分数据集
sccForestData<-rbind(sstselfEval("trainPo","Training"),
                     sstselfEval("testPo","Test"),sstselfEval("independent","Independent"))
rownames(sccForestData)<-paste(c(rep("Training",5),rep("Test",5),rep("Independent",5)),
                               rownames(sccForestData),sep="-")
sccForestData<-as.data.frame(sccForestData)
Data<-data.frame(lnHR=sccForestData$lnHR,
                 selnHR=sccForestData$selnHR)
rownames(Data)<-rownames(sccForestData)
pdf("sccForest.pdf",width=6,height=4)
forest(metagen(Data$lnHR,Data$selnHR,sm="HR"),overall = F,squaresize = 0.4,
       weight.study = "same",xlim = c(0.5,4),leftcols = "studlab",col.square="#2D5662",
       studlab=rownames(Data),rightcols = c("effect.ci"),ref=1,plotwidth = "3cm",
       col.predict.lines="#2D5662")
dev.off()
HR<-sccForestData$HR
lnHR<-sccForestData$lnHR
selnHR<-sccForestData$selnHR
low.HR<-exp(lnHR-1.96*selnHR)
high.HR<-exp(lnHR+1.96*selnHR)
p<-sccForestData$`p value`
a<-data.frame(sccForestData$`Low Risk`,
              sccForestData$`High Risk`,
              format(HR,digits=2),
              paste(format(low.HR,digits=2),
                    format(high.HR,digits=2),sep="-"),
              p=format(p,digits=2,scientific = 2))
colnames(a)<-c("Low Risk","High Risk","HR","95%CI of HR","P value")
rownames(a)<-rownames(sccForestData)
write.csv(a,"sccForestData.csv")

#global
ooo<-sstselfEval("global","global")
sccForestData<-ooo[-1,]
sccForestData<-as.data.frame(sccForestData)
HR<-sccForestData$HR
lnHR<-sccForestData$lnHR
selnHR<-sccForestData$selnHR
low.HR<-exp(lnHR-1.96*selnHR)
high.HR<-exp(lnHR+1.96*selnHR)
p<-sccForestData$`p value`
a<-data.frame(sccForestData$`Low Risk`,
              sccForestData$`High Risk`,
              round(HR,digits=2),
              paste(round(low.HR,digits=2),
                    round(high.HR,digits=2),sep="-"),
              p=format(p,digits=2,scientific = 2))
colnames(a)<-c("Low Risk","High Risk","HR","95%CI of HR","P value")
rownames(a)<-rownames(sccForestData)
write.csv(a,"sccForestData.csv")

#分训练和测试
sccForestData<-rbind(sstselfEval("trainPo","Training"),
                     sstselfEval("TestALL","TestALL"))
######离散分数单多cox#####
unicox<-rbind(sstselfEvalUM("trainPo")[[1]],
              sstselfEvalUM("testPo")[[1]],
              sstselfEvalUM("independent")[[1]])
rownames(unicox)<-paste(c(rep("Training",3),rep("Test",3),rep("Independent",3)),
                               rownames(unicox),sep="-")
write.csv(createTable(as.data.frame(unicox)),"selfUnicox.csv",row.names = F)
multicox<-rbind(sstselfEvalUM("trainPo")[[2]],
              sstselfEvalUM("testPo")[[2]],
              sstselfEvalUM("independent")[[2]])
rownames(multicox)<-paste(c(rep("Training",3),rep("Test",3),rep("Independent",3)),
                          rownames(multicox),sep="-")
write.csv(createTable(as.data.frame(multicox)),"selfmulticox.csv",row.names = F)
multicox<-as.data.frame(multicox)
Data<-data.frame(lnHR=multicox$coef,
                 selnHR=multicox$`se(coef)`)
rownames(Data)<-rownames(multicox)
pdf("selfMulticoxForest.pdf",width=6,height=4)
forest(metagen(Data$lnHR,Data$selnHR,sm="HR"),overall = F,squaresize = 0.4,
       weight.study = "same",xlim = c(0.5,2.5),leftcols = "studlab",col.square="#2D5662",
       studlab=rownames(Data),rightcols = c("effect.ci"),ref=1,plotwidth = "3cm",
       col.predict.lines="#2D5662")
dev.off()
######连续分数单多cox######
#扔掉
unicox<-rbind(sstselfEvalUMContinue("trainPo")[[1]],
              sstselfEvalUMContinue("testPo")[[1]],
              sstselfEvalUMContinue("independent")[[1]])
rownames(unicox)<-paste(c(rep("Training",3),rep("Test",3),rep("Independent",3)),
                        rownames(unicox),sep="-")
#####6.单多cox自我评估#####
gradeStatus<-c(0,rep(1,5),0,rep(1,8),2)
#获取训练集各个数据集的stage,grade,age状态
clinicalTrain<-clinical[trainPo]
gradeStatusTrain<-gradeStatus[trainPo]
result<-HGClinicalStatus(clinicalTrain,gradeStatusTrain)
idx.tumorstage.train<-result[,1]
idx.tumorgrade.train<-result[,2]
idx.age.train<-result[,3]
stageNameTrain<-c(rep("tumorstage",15),"tumour_stage")[trainPo]
gradeNameTrain<-c(rep("grade",15),"tumour_grade")[trainPo]
ageNameTrain<-c(rep("age_at_initial_pathologic_diagnosis",15),
           "donor_age_at_diagnosis")[trainPo]
timeNameTrain<-c(rep("days_to_death",15),"donor_survival_time")[trainPo]
eventNameTrain<-c(rep("vital_status",15),"donor_vital_status")[trainPo]
selfEvalResult<-list()
par(mfrow=c(4,4))#要手动存的哟,名字：model ROC Train
for(i in 1:length(clinicalTrain)){#length(clinicalTrain)
  score.i<-scale(as.numeric(ssgsea.out[[i]][4,]))
  phe<-clinicalTrain[[i]]
  
  sCol<-which(colnames(phe)==stageNameTrain[i])
  gCol<-which(colnames(phe)==gradeNameTrain[i])
  aCol<-which(colnames(phe)==ageNameTrain[i])
  timeCol<-which(colnames(phe)==timeNameTrain[i])
  eventCol<-which(colnames(phe)==eventNameTrain[i])
  
  result<-selfEval(score.i,phe,names(clinicalTrain)[i],
                   idx.tumorstage.train[i],idx.tumorgrade.train[i],
                   idx.age.train[i],sCol,gCol,aCol,timeCol,eventCol)
  selfEvalResult[[names(clinicalTrain)[i]]]<-result
}
#cox结果
selfEvalCoxResult<-c()
for(i in 1:length(selfEvalResult)){
  x<-selfEvalResult[[i]]
  Xname<-names(selfEvalResult)[i]
  selfEvalCoxResult<-rbind(selfEvalCoxResult,combineCoxResult(x,Xname))
}
write.csv(selfEvalCoxResult,"selfEvalCoxResultTrain.csv",row.names = F)
#roc.ci
threeYearSurvival<-sapply(selfEvalResult,function(x){x[[3]][1]})
fiveYearSurvival<-sapply(selfEvalResult,function(x){x[[3]][2]})
c.index<-sapply(selfEvalResult,function(x){x[[3]][3]})
sd.c.index<-sapply(selfEvalResult,function(x){x[[3]][4]})
low.ci<-c.index-1.96*sd.c.index
high.ci<-c.index+1.96*sd.c.index

roc.ci<-data.frame(names(selfEvalResult),format(threeYearSurvival,digits = 4),
                   format(fiveYearSurvival,digits = 4),
                   format(c.index,digits = 4),
                   paste(format(low.ci,digits = 4),
                         format(high.ci,digits = 4),sep="-"))
names(roc.ci)<-c("Datasets","Three-year survival",
                 "Five-year survival","C-index","95%CI of C-index")
write.csv(roc.ci,"roc_ciTrain.csv",row.names = F)
#forest.multi
pdf("multiForestTrain.pdf",width = 7,height = 7)
for(i in 1:length(selfEvalResult)){
  a<-selfEvalResult[[i]][[1]]
  rownames(a)[grepl("sCo",rownames(a),fixed = T)]<-"stage"
  rownames(a)[grepl("aCol",rownames(a),fixed = T)]<-"age"
  rownames(a)[grepl("gCol",rownames(a),fixed = T)]<-"grade"
  a[,5]<-round(a[,5],digits = 4)
  colnames(a)<-c("lnhr","hr","se_lnhr","z","p")
  a<-as.data.frame(a)
  meta.i<-metagen(lnhr,se_lnhr,sm="HR",data = a,studlab = rownames(a))
  forest(meta.i, leftcols = c("studlab"),leftlabs = c("Study Element"),
         rightcols = c("effect.ci","p"),fontsize = 13,
         rightlabs = "p",just.studlab = "center",
         overall = F,col.square = "#363886",lwd=2,
         weight.study = "same",col.inside="black",xlab=names(selfEvalResult)[i])
}
dev.off()
#forest.uni
HR<-sapply(selfEvalResult,function(x){x[[2]][1,2]})
lnHR<-sapply(selfEvalResult,function(x){x[[2]][1,1]})
selnHR<-sapply(selfEvalResult,function(x){x[[2]][1,3]})
low.HR<-exp(lnHR-1.96*selnHR)
high.HR<-exp(lnHR+1.96*selnHR)
C_index<-sapply(selfEvalResult,function(x){x[[3]][3]})
se_C_index<-sapply(selfEvalResult,function(x){x[[3]][4]})
low_c_index<-C_index-1.96*se_C_index
high_c_index<-C_index+1.96*se_C_index
HR_C<-data.frame(HR,lnHR,selnHR,low.HR,high.HR,
                 C_index,se_C_index,low_c_index,high_c_index)
write.csv(HR_C,"HR_CTrain.csv")

lnHR<-meta.summaries(HR_C[,2],HR_C[,3],method="fixed",logscale=T)
if(lnHR$het[3]>0.05){
  summary_HR<-c(lnHR$summary,lnHR$se.summary)
}else{
  lnHR<-meta.summaries(HR_C[,2],HR_C[,3],method="random",logscale=T)
  summary_HR<-c(lnHR$summary,lnHR$se.summary)
}
#求C-statics的meta结果
C<-meta.summaries(HR_C[,6],HR_C[,7],method="fixed",logscale=F)
if(C$het[3]>0.05){
  summary_C<-c(C$summary,C$se.summary)
}else{
  C<-meta.summaries(HR_C[,6],HR_C[,7],method="random",logscale=F)
  summary_C<-c(C$summary,C$se.summary)
}
ssgsea_colla_Train<-c(summary_HR,summary_C)
write.csv(ssgsea_colla_Train,"ssgsea_colla_Train.csv")
ol<-colorRampPalette(c("#363886","#9ECAD4"))
pdf("unicox-HR(no ln) Train.pdf ",width=6,height=5)
a<-metagen(HR_C$lnHR,HR_C$selnHR,sm="HR")
data<-data.frame(HR=c(HR_C$lnHR,a$TE.random),
                 seHR=c(HR_C$selnHR,a$seTE.random))
rownames(data)<-c(rownames(HR_C),"Total")
forest(metagen(data$HR,data$seHR,sm="HR"),overall = F,squaresize = 0.4,
       weight.study = "same",xlim = c(1,2.5),leftcols = "studlab",col.square="#363886",
       studlab=rownames(data),rightcols = c("effect.ci"),ref=1,plotwidth = "3cm",
       at=c(1,1.5,2,2.5))
dev.off()
pdf("uniCox-HR Train.pdf",width = 10,height = 7)
a<-metagen(HR_C$lnHR,HR_C$selnHR,sm="lnHR")
data<-data.frame(lnHR=c(HR_C$lnHR,a$TE.random),
                 selnHR=c(HR_C$selnHR,a$seTE.random))
rownames(data)<-c(rownames(HR_C),"Total")
forest(metagen(data$lnHR,data$selnHR,sm="lnHR"),overall = F,squaresize = 0.4,
       weight.study = "same",xlim = c(0,1),leftcols = "studlab",col.square="#363886",
      studlab=rownames(data),rightcols = c("effect.ci"),at=c(0,0.2,0.4,0.6,0.8,1),plotwidth="3.5cm")
dev.off()

pdf("uniCox-C-index Train.pdf",width = 10,height = 7)
a<-metagen(HR_C$C_index,HR_C$se_C_index,sm="C-index")
data<-data.frame(C_index=c(HR_C$C_index,a$TE.random),
                 se_C_index=c(HR_C$se_C_index,a$seTE.random))
rownames(data)<-c(rownames(HR_C),"Total")
forest(metagen(data$C_index,data$se_C_index,sm="C-index"),overall = F,squaresize = 0.4,
       weight.study = "same",xlim = c(0.5,0.8),leftcols = "studlab",col.square="#363886",
       studlab=rownames(data),rightcols = c("effect.ci"),ref=0.5,plotwidth = "3cm",
       at=c(0.5,0.6,0.7,0.8))
dev.off()
#####7.单多cox测单评估#####
gradeStatus<-c(0,rep(1,5),0,rep(1,8),2)
#获取训练集各个数据集的stage,grade,age状态
clinicalTest<-clinical[-trainPo]
gradeStatusTest<-gradeStatus[-trainPo]
result<-HGClinicalStatus(clinicalTest,gradeStatusTest)
idx.tumorstage.test<-result[,1]
idx.tumorgrade.test<-result[,2]
idx.age.test<-result[,3]
stageNameTest<-c(rep("tumorstage",15),"tumour_stage")[-trainPo]
gradeNameTest<-c(rep("grade",15),"tumour_grade")[-trainPo]
ageNameTest<-c(rep("age_at_initial_pathologic_diagnosis",15),
                "donor_age_at_diagnosis")[-trainPo]
timeNameTest<-c(rep("days_to_death",15),"donor_survival_time")[-trainPo]
eventNameTest<-c(rep("vital_status",15),"donor_vital_status")[-trainPo]
selfEvalResult<-list()
par(mfrow=c(2,2))#要手动存的哟,名字：model ROC Test
for(i in 1:length(clinicalTest)){
  score.i<-scale(as.numeric(ssgsea.out.test[[i]][4,]))
  phe<-clinicalTest[[i]]
  
  sCol<-which(colnames(phe)==stageNameTest[i])
  gCol<-which(colnames(phe)==gradeNameTest[i])
  aCol<-which(colnames(phe)==ageNameTest[i])
  timeCol<-which(colnames(phe)==timeNameTest[i])
  eventCol<-which(colnames(phe)==eventNameTest[i])
  
  result<-selfEval(score.i,phe,names(clinicalTest)[i],
                   idx.tumorstage.test[i],idx.tumorgrade.test[i],
                   idx.age.test[i],sCol,gCol,aCol,timeCol,eventCol)
  selfEvalResult[[names(clinicalTest)[i]]]<-result
}
#cox结果
selfEvalCoxResult<-c()
for(i in 1:length(selfEvalResult)){
  x<-selfEvalResult[[i]]
  Xname<-names(selfEvalResult)[i]
  selfEvalCoxResult<-rbind(selfEvalCoxResult,combineCoxResult(x,Xname))
}
write.csv(selfEvalCoxResult,"selfEvalCoxResultTest.csv",row.names = F)
#roc.ci
threeYearSurvival<-sapply(selfEvalResult,function(x){x[[3]][1]})
fiveYearSurvival<-sapply(selfEvalResult,function(x){x[[3]][2]})
c.index<-sapply(selfEvalResult,function(x){x[[3]][3]})
sd.c.index<-sapply(selfEvalResult,function(x){x[[3]][4]})
low.ci<-c.index-1.96*sd.c.index
high.ci<-c.index+1.96*sd.c.index

roc.ci<-data.frame(names(selfEvalResult),format(threeYearSurvival,digits = 4),
                   format(fiveYearSurvival,digits = 4),
                   format(c.index,digits = 4),
                   paste(format(low.ci,digits = 4),
                         format(high.ci,digits = 4),sep="-"))
names(roc.ci)<-c("Datasets","Three-year survival",
                 "Five-year survival","C-index","95%CI of C-index")
write.csv(roc.ci,"roc_ci_Test.csv",row.names = F)
#forest.multi
pdf("multiForestTest.pdf",width = 7,height = 7)
for(i in 1:length(selfEvalResult)){
  a<-selfEvalResult[[i]][[1]]
  rownames(a)[grepl("sCo",rownames(a),fixed = T)]<-"stage"
  rownames(a)[grepl("aCol",rownames(a),fixed = T)]<-"age"
  rownames(a)[grepl("gCol",rownames(a),fixed = T)]<-"grade"
  a[,5]<-round(a[,5],digits = 4)
  colnames(a)<-c("lnhr","hr","se_lnhr","z","p")
  a<-as.data.frame(a)
  meta.i<-metagen(lnhr,se_lnhr,sm="HR",data = a,studlab = rownames(a))
  forest(meta.i, leftcols = c("studlab"),leftlabs = c("Study Element"),
         rightcols = c("effect.ci","p"),fontsize = 13,
         rightlabs = "p",just.studlab = "center",
         overall = F,col.square = "#363886",lwd=2,
         weight.study = "same",col.inside="black",xlab=names(selfEvalResult)[i])
}
dev.off()
#forest.uni
HR<-sapply(selfEvalResult,function(x){x[[2]][1,2]})
lnHR<-sapply(selfEvalResult,function(x){x[[2]][1,1]})
selnHR<-sapply(selfEvalResult,function(x){x[[2]][1,3]})
low.HR<-exp(lnHR-1.96*selnHR)
high.HR<-exp(lnHR+1.96*selnHR)
C_index<-sapply(selfEvalResult,function(x){x[[3]][3]})
se_C_index<-sapply(selfEvalResult,function(x){x[[3]][4]})
low_c_index<-C_index-1.96*se_C_index
high_c_index<-C_index+1.96*se_C_index
HR_C<-data.frame(HR,lnHR,selnHR,low.HR,high.HR,
                 C_index,se_C_index,low_c_index,high_c_index)
write.csv(HR_C,"HR_C_Test.csv")

lnHR<-meta.summaries(HR_C[,2],HR_C[,3],method="fixed",logscale=T)
if(lnHR$het[3]>0.05){
  summary_HR<-c(lnHR$summary,lnHR$se.summary)
}else{
  lnHR<-meta.summaries(HR_C[,2],HR_C[,3],method="random",logscale=T)
  summary_HR<-c(lnHR$summary,lnHR$se.summary)
}
#求C-statics的meta结果
C<-meta.summaries(HR_C[,6],HR_C[,7],method="fixed",logscale=F)
if(C$het[3]>0.05){
  summary_C<-c(C$summary,C$se.summary)
}else{
  C<-meta.summaries(HR_C[,6],HR_C[,7],method="random",logscale=F)
  summary_C<-c(C$summary,C$se.summary)
}
ssgsea_colla_Test<-c(summary_HR,summary_C)
write.csv(ssgsea_colla_Test,"ssgsea_colla_Test.csv")
ol<-colorRampPalette(c("#363886","#9ECAD4"))

pdf("uniCox-HR Test.pdf",width = 10,height = 7)
forest(metagen(HR_C$lnHR,HR_C$selnHR,sm="HR"),
       col.square =  ol(nrow(HR_C)),col.diamond = "#70C1B3",
       studlab=rownames(HR_C))
dev.off()

pdf("uniCox-C-index Test .pdf",width = 10,height = 7)
forest(metagen(HR_C$C_index,HR_C$se_C_index,sm="C"),
       xlim = c(0,1),col.square =  ol(nrow(HR_C)),col.diamond = "#70C1B3",
       studlab=rownames(HR_C),ref = 0.5)
dev.off()
#####8.训练集整体+测试集整体评估#####
clinicalTrain<-clinical[trainPo]
clinicalTest<-clinical[-trainPo]
clinicalMatrixTrain<-getSurvivalMatrix(clinicalTrain,"global","Patch et al.")
clinicalMatrixTest<-getSurvivalMatrix(clinicalTest,"global","Patch et al.")
par(mfrow=c(3,4))#name:globalROC
globalTrain<-globalEval(scale(score),"global Training",clinicalMatrixTrain,1,2)
globalTest<-globalEval(scale(score.test),"global Test",clinicalMatrixTest,1,2)
global<-rbind(globalTrain,globalTest)
#####8.1训整+测整1+测整2评估#####
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

pdf("global Trainingsst.pdf",width=6,height=6)
globalTrain<-globalEval(scale(score),"global Trainingsst",clinicalMatrixTrain,1,2)
dev.off()
pdf("global Testsst.pdf",width=6,height=6)
globalTest<-globalEval(scale(score.test),"global Testsst",clinicalMatrixTest,1,2)
dev.off()
pdf("global independsst.pdf",width=6,height=6)
globalIndepend<-globalEval(scale(score.test.independ),
                           "global independsst",clinicalMatrixIndepend,1,2)
dev.off()
#####所有样本的生存ROC
clinicalAll<-rbind(clinicalMatrixTrain,clinicalMatrixTest,clinicalMatrixIndepend)
score.all<-c(score,score.test,score.test.independ)
pdf("global all samplesst.pdf",width=6,height=6)
globalall<-globalEval(scale(score.all),"global all samplesst",clinicalAll,1,2)
dev.off()
######折线图####
data<-as.data.frame(global[,7:11])
data1<-data.frame(value=unlist(data),
                  survivalTime=rep(1:5,each=2),
                  class=rep(rownames(data),times=5))
ggplot(data1,aes(survivalTime,value,group=class))+
  geom_line(size=1,aes(color=class))+
  geom_point(size=2.5,aes(color=class))+
  geom_hline(yintercept = 0.5,linetype="dashed")+
  theme_bw()+
  scale_color_manual(values=c("#E07A9A","#805BA6"),
                     labels=c("Test set","Training set"))+
  labs(x="Survival time (Years)",y="AUC",title="Time depedent ROC curve")+
  theme(panel.border = element_blank(),panel.grid=element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(color = "black",size=15),
        axis.text.y=element_text(color="black",size=15),
        plot.title = element_text(colour = "black", 
                                  size = 15, vjust = 1,hjust = 0.5),
        axis.title = element_text(size=15),
        legend.title = element_blank(),
        legend.justification=c(0,0),
        legend.position=c(0.01,0.05),
        legend.text=element_text(size=12,face="bold"))
ggsave("global 折线.pdf",width=6,height=6)
##接折线图上方程序##
global<-as.data.frame(global)
global_Value_se<-data.frame(global$lnHR,global$selnHR,
                            global$C_index,global$se_C_index)
write.csv(global_Value_se,"global_Value_se.csv")
global1<-data.frame(rownames(global),
                    format(global$HR,digits=4),
                    paste(format(global$low.HR,digits=4),
                          format(global$high.HR,digits=4),sep="-"),
                    format(global$p.HR,digits=3,scientific = T),
                    format(global$threeYearSurvival,digits=4),
                    format(global$fiveYearSurvival,digits=4),
                    format(global$C_index,digits=4),
                    paste(format(global$low_C_index,digits=4),
                          format(global$high_C_index,digits=4),sep="-"))
names(global1)<-c("Variables","HR","95%CI of HR","P value of HR",
                  "three-Year Survival",
                  "five-Year Survival","C-index","95%CI of C-index")
write.csv(global1,"global.csv",row.names = F)

pdf("global HR Forest.pdf",width = 7,height = 7)
global$p.HR<-format(global$p.HR,digits = 3,scientific = T)
meta.i<-metagen(lnHR,selnHR,sm="HR",data = global,studlab = rownames(global))
forest(meta.i, leftcols = c("studlab"),leftlabs = c("Study Element"),
       rightcols = c("effect.ci","p.HR"),fontsize = 13,just.studlab = "center",
       overall = F,col.square = "#363886",lwd=2,rightlabs = "P value",
       weight.study = "same",col.inside="black")  
dev.off()

pdf("global C-index Forest.pdf",width = 7,height = 7)
meta.i<-metagen(C_index,se_C_index,sm="C-index",data = global,
                studlab = rownames(global))
forest(meta.i, leftcols = c("studlab"),leftlabs = c("Study Element"),
       rightcols = c("effect.ci"),fontsize = 13,just.studlab = "center",
       overall = F,col.square = "#363886",lwd=2,xlim = c(0,1),
       weight.study = "same",col.inside="black",ref = 0.5)  

dev.off()

####9.与其他模型比较####
####9.1 14个模型统计量获取####
#模型相似度

source('G:/haodapeng-data/deleteScatter/code/lower2symmetric.R')
source('G:/haodapeng-data/deleteScatter/codes/signatureSimility.R')
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
write.csv(ModelGeneSim,"ModelGeneSim.csv")
write.csv(ModelGeneSimP,"ModelGeneSimP.csv")
#训单
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
save(predictive_ability_Train,file="predictive_ability_Train.Rdata")
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
write.csv(predictive_ability_Train_excel,"14 Model seperate Performance Train.csv")
#训单合并
load("predictive_ability_Train.Rdata")
matrix_meta_Train<-c()
for(i in 1:nrow(predictive_ability_Train[[1]])){
  matrix_premeta<-c()
  for(j in 1:length(predictive_ability_Train)){
    matrix_premeta<-rbind(matrix_premeta,predictive_ability_Train[[j]][i,])
  }
  #求HR的meta结果。
  lnHR<-meta.summaries(matrix_premeta[,1],matrix_premeta[,2],method="fixed",logscale=T)
  if(lnHR$het[3]>0.05){
    summary_HR<-c(lnHR$summary,lnHR$se.summary)
  }else{
    lnHR<-meta.summaries(matrix_premeta[,1],matrix_premeta[,2],method="random",logscale=T)
    summary_HR<-c(lnHR$summary,lnHR$se.summary)
  }
  #求C-statics的meta结果
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
write.csv(matrix_meta_Train,"matrix_meta_Train.csv")
#训整
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
write.csv(ForteenPerfoTrain,"ForteenPerfoTrain.csv")

#测单
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
write.csv(predictive_ability_Test_excel,"14 Model seperate Performance Test.csv")
#测单合并
load("predictive_ability_Test.Rdata")
matrix_meta_Test<-c()
for(i in 1:nrow(predictive_ability_Test[[1]])){
  matrix_premeta<-c()
  for(j in 1:length(predictive_ability_Test)){
    matrix_premeta<-rbind(matrix_premeta,predictive_ability_Test[[j]][i,])
  }
  #求HR的meta结果。
  lnHR<-meta.summaries(matrix_premeta[,1],matrix_premeta[,2],method="fixed",logscale=T)
  if(lnHR$het[3]>0.05){
    summary_HR<-c(lnHR$summary,lnHR$se.summary)
  }else{
    lnHR<-meta.summaries(matrix_premeta[,1],matrix_premeta[,2],method="random",logscale=T)
    summary_HR<-c(lnHR$summary,lnHR$se.summary)
  }
  #求C-statics的meta结果
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
write.csv(matrix_meta_Test,"matrix_meta_Test.csv")

#测整
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
write.csv(ForteenPerfoTest,"ForteenPerfoTest.csv")
####9.2 所有模型(15个)比较####
####9.2.1训单####
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
####9.2.2训单合并####
ssgsea_colla_Train<-as.numeric(unlist(read.csv("ssgsea_colla_Train.csv",
                                               header=T,row.names = 1,as.is = T)))
matrix_meta_Train<-read.csv("matrix_meta_Train.csv",header=T,row.names = 1,as.is = T)#训单合并
colnames(matrix_meta_Train)<-names(ssgsea_colla_Train)<-c("lnHR","selnHR",
                                                    "C_index","se_C_index")
                                                                                      
Data<-rbind(matrix_meta_Train,ssgsea_colla_Train)
rownames(Data)<-c(rownames(matrix_meta_Train),"-FIPI")

pdf("14 VS ssgsea CombineSeperate trainC-index.pdf",width=7,height = 7)
Data<-as.data.frame(Data[order(Data[,3],decreasing = T),])
meta.i<-metagen(C_index,se_C_index,sm="C-index",data = Data,
                studlab = rownames(Data))
forest(meta.i, leftcols = c("studlab"),leftlabs = c("Study Element"),
       rightcols = c("effect.ci"),fontsize = 13,just.studlab = "center",
       overall = F,col.square = "#363886",lwd=2,xlim = c(0.45,0.7),plotwidth="3.5cm",
       weight.study = "same",col.inside="black",ref = 0.5,squaresize = 0.5)
dev.off()

pdf("14 VS ssgsea CombineSeperate trainHR.pdf",width=7,height = 7)
Data<-as.data.frame(Data[order(Data[,1],decreasing = T),])
meta.i<-metagen(lnHR,selnHR,sm="lnHR",data = Data,
                studlab = rownames(Data))
forest(meta.i, leftcols = c("studlab"),leftlabs = c("Study Element"),
       rightcols = c("effect.ci"),fontsize = 13,just.studlab = "center",
       overall = F,col.square = "#363886",lwd=2,plotwidth="3.5cm",squaresize = 0.5,
       weight.study = "same",col.inside="black",ref = 0,xlim = c(-0.2,0.6))
dev.off()
####9.2.3训整####
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

pdf("14 VS ssgsea globalTrain HR.pdf",width=7,height = 7)
Data<-as.data.frame(Data[order(Data[,1],decreasing = T),])
meta.i<-metagen(lnHR,selnHR,sm="HR",data = Data,
                studlab = rownames(Data))
forest(meta.i, leftcols = c("studlab"),leftlabs = c("Study Element"),
       rightcols = c("effect.ci"),fontsize = 13,just.studlab = "center",
       overall = F,col.square = "#363886",lwd=2,xlim = c(0.9,7),
       weight.study = "same",col.inside="black",ref = 1)
dev.off()
####9.2.4测单####
HR_C_Test<-read.csv("HR_C_Test.csv",header=T,row.names = 1,as.is = T)
load("predictive_ability_Test.Rdata")
pdf("14 VS ssgsea seperate test.pdf",width = 7,height = 7)
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
dev.off()
####9.2.5测合并####
ssgsea_colla_Test<-as.numeric(unlist(read.csv("ssgsea_colla_Test.csv",
                                               header=T,row.names = 1,as.is = T)))
matrix_meta_Test<-read.csv("matrix_meta_Test.csv",header=T,row.names = 1,as.is = T)#训单合并
colnames(matrix_meta_Test)<-names(ssgsea_colla_Test)<-c("lnHR","selnHR",
                                                          "C_index","se_C_index")

Data<-rbind(matrix_meta_Test,ssgsea_colla_Test)
rownames(Data)<-c(rownames(ForteenPerfoTest),"ssgsea")

pdf("14 VS ssgsea CombineSeperate testC-index.pdf",width=7,height = 7)
Data<-as.data.frame(Data[order(Data[,3],decreasing = T),])
meta.i<-metagen(C_index,se_C_index,sm="C-index",data = Data,
                studlab = rownames(Data))
forest(meta.i, leftcols = c("studlab"),leftlabs = c("Study Element"),
       rightcols = c("effect.ci"),fontsize = 13,just.studlab = "center",
       overall = F,col.square = "#363886",lwd=2,xlim = c(0.4,0.7),
       weight.study = "same",col.inside="black",ref = 0.5)
dev.off()

pdf("14 VS ssgsea CombineSeperate testHR.pdf",width=7,height = 7)
Data<-as.data.frame(Data[order(Data[,1],decreasing = T),])
meta.i<-metagen(lnHR,selnHR,sm="HR",data = Data,
                studlab = rownames(Data))
forest(meta.i, leftcols = c("studlab"),leftlabs = c("Study Element"),
       rightcols = c("effect.ci"),fontsize = 13,just.studlab = "center",
       overall = F,col.square = "#363886",lwd=2,
       weight.study = "same",col.inside="black",ref = 1)
dev.off()
#####9.2.6测整####
colnames(global_Value_se)<-colnames(ForteenPerfoTest)<-c("lnHR","selnHR",
                                                          "C_index","se_C_index")
Data<-rbind(ForteenPerfoTest,global_Value_se[2,])
rownames(Data)<-c(rownames(ForteenPerfoTest),"ssgsea")
pdf("14 VS ssgsea globalTestC-index.pdf",width=7,height = 7)
Data<-as.data.frame(Data[order(Data[,3],decreasing = T),])
meta.i<-metagen(C_index,se_C_index,sm="C-index",data = Data,
                studlab = rownames(Data))
forest(meta.i, leftcols = c("studlab"),leftlabs = c("Study Element"),
       rightcols = c("effect.ci"),fontsize = 13,just.studlab = "center",
       overall = F,col.square = "#363886",lwd=2,xlim = c(0.4,0.7),
       weight.study = "same",col.inside="black",ref = 0.5)
dev.off()

pdf("14 VS ssgsea globalTest HR.pdf",width=7,height = 7)
Data<-as.data.frame(Data[order(Data[,1],decreasing = T),])
meta.i<-metagen(lnHR,selnHR,sm="HR",data = Data,
                studlab = rownames(Data))
forest(meta.i, leftcols = c("studlab"),leftlabs = c("Study Element"),
       rightcols = c("effect.ci"),fontsize = 13,just.studlab = "center",
       overall = F,col.square = "#363886",lwd=2,xlim = c(0.9,3),
       weight.study = "same",col.inside="black",ref = 1)
dev.off()

#####9.2.7 模型jaccard箱式图#####
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
write.csv(jaccardMatrix_model,"jaccardMatrix_model.csv")

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
ggsave("model sim based on jaccard.pdf",width=6,height=4)
#####10.所有数据集的hebing比较#####
#####10.1 ssgsea#####
gradeStatus<-c(0,rep(1,5),0,rep(1,8),2)
gradeStatus[14]<-2
#获取训练集各个数据集的stage,grade,age状态
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
#par(mfrow=c(4,4))#要手动存的哟,名字：model ROC Train
for(i in 1:length(clinical)){#单、多cox都是基于连续分值做的
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
####unicox
uni_ssgsea_colla<-unicoxForestData(selfEvalResult)[[1]]
####multicox
multi_ssgsea_colla<-MulticoxForestData(selfEvalResult)[[1]]

#####10.2 others#####
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
  #求HR的meta结果。
  lnHR<-meta.summaries(matrix_premeta[,1],matrix_premeta[,2],method="fixed",logscale=T)
  if(lnHR$het[3]>0.05){
    summary_HR<-c(lnHR$summary,lnHR$se.summary)
  }else{
    lnHR<-meta.summaries(matrix_premeta[,1],matrix_premeta[,2],method="random",logscale=T)
    summary_HR<-c(lnHR$summary,lnHR$se.summary)
  }
  #求C-statics的meta结果
  C<-meta.summaries(matrix_premeta[,3],matrix_premeta[,4],method="fixed",logscale=F)
  if(C$het[3]>0.05){
    summary_C<-c(C$summary,C$se.summary)
  }else{
    C<-meta.summaries(matrix_premeta[,3],matrix_premeta[,4],method="random",logscale=F)
    summary_C<-c(C$summary,C$se.summary)
  }
  uni_matrix_meta<-rbind(uni_matrix_meta,c(summary_HR,summary_C))
}
rownames(uni_matrix_meta)<-rownames(uni_predictive_ability[[1]])
colnames(uni_matrix_meta)<-c("lnHR","se(lnHR)","C-index","se(C-index)")
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
  #求HR的meta结果。
  lnHR<-meta.summaries(matrix_premeta[,1],matrix_premeta[,2],method="fixed",logscale=T)
  if(lnHR$het[3]>0.05){
    summary_HR<-c(lnHR$summary,lnHR$se.summary)
  }else{
    lnHR<-meta.summaries(matrix_premeta[,1],matrix_premeta[,2],method="random",logscale=T)
    summary_HR<-c(lnHR$summary,lnHR$se.summary)
  }
  #求C-statics的meta结果
  C<-meta.summaries(matrix_premeta[,3],matrix_premeta[,4],method="fixed",logscale=F)
  if(C$het[3]>0.05){
    summary_C<-c(C$summary,C$se.summary)
  }else{
    C<-meta.summaries(matrix_premeta[,3],matrix_premeta[,4],method="random",logscale=F)
    summary_C<-c(C$summary,C$se.summary)
  }
  Multi_matrix_meta<-rbind(Multi_matrix_meta,c(summary_HR,summary_C))
}
rownames(Multi_matrix_meta)<-rownames(multi_predictive_ability[[1]])
colnames(Multi_matrix_meta)<-c("lnHR","se(lnHR)","C-index","se(C-index)")
#####10.3 forest#####
#####uni####
Data<-rbind(uni_matrix_meta,uni_ssgsea_colla)
rownames(Data)<-c(rownames(uni_matrix_meta),"ssgsea")
rownames(Data)[which("ours116"==rownames(Data))]<-"lu14"
colnames(Data)<-c("lnHR","selnHR","C.index","se.C.index")

pdf("uni 14 VS ssgsea CombineSeperate trainHR.pdf",width=7,height = 7)
Data<-as.data.frame(Data[order(Data[,1],decreasing = T),])
uniHRname<-rownames(Data)
meta.i<-metagen(lnHR,selnHR,sm="HR",data = Data,
                studlab = rownames(Data))
forest(meta.i, leftcols = c("studlab"),leftlabs = c("Study Element"),
       rightcols = c("effect.ci"),fontsize = 13,just.studlab = "center",
       at=c(0.9,1,1.1,1.2,1.3,1.4,1.5,1.6),
       overall = F,col.square = "#363886",lwd=2,squaresize = 0.4,xlim=c(0.9,1.6),
       weight.study = "same",col.inside="black",ref = 1,plotwidth="3.5cm")
dev.off()

pdf("uni 14 VS ssgsea CombineSeperate trainC-index.pdf",width=7,height = 7)
Data<-as.data.frame(Data[order(Data[,3],decreasing = T),])
uniCname<-rownames(Data)
meta.i<-metagen(C.index,se.C.index,sm="C-index",data = Data,
                studlab = rownames(Data))
forest(meta.i, leftcols = c("studlab"),leftlabs = c("Study Element"),
       rightcols = c("effect.ci"),fontsize = 13,just.studlab = "center",
       overall = F,col.square = "#363886",lwd=2,xlim = c(0.4,0.7),
       at=c(0.4,0.5,0.6,0.7),plotwidth="3.5cm",squaresize = 0.4,
       weight.study = "same",col.inside="black",ref = 0.5)
dev.off()


Data<-Data[c(15,3,1,6,5,2,7),]
Data<-as.data.frame(Data[order(Data[,3],decreasing = T),])
Data$model<-factor(rownames(Data),levels =rownames(Data) )
ggplot(Data, aes(x=model, y=C.index)) + 
  geom_bar(stat ="identity",width = 0.8,fill="#BE0005") +#"#EFC000"
  geom_errorbar(aes(ymin=C.index-1.96*se.C.index, ymax=C.index+1.96*se.C.index), width=.1) +
  #geom_errorbar(aes(ymin=aucOf3-se.aucof3, ymax=aucOf3+se.aucof3), width=.1)
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
ggsave("coxmodels_average_C.pdf",width = 7,height=7)

#####multi#####
Data<-rbind(Multi_matrix_meta,multi_ssgsea_colla)
rownames(Data)<-c(rownames(Multi_matrix_meta),"-FIPI")
rownames(Data)[which("ours116"==rownames(Data))]<-"lu14"
colnames(Data)<-c("lnHR","selnHR","C.index","se.C.index")

pdf("multi 14 VS ssgsea CombineSeperate trainHR.pdf",width=7,height = 7)
Data<-as.data.frame(Data[order(Data[,1],decreasing = T),])
MultiHRname<-rownames(Data)
meta.i<-metagen(lnHR,selnHR,sm="HR",data = Data,
                studlab = rownames(Data))
forest(meta.i, leftcols = c("studlab"),leftlabs = c("Study Element"),
       rightcols = c("effect.ci"),fontsize = 13,just.studlab = "center",
       at=c(0.9,1,1.1,1.2,1.3,1.4,1.6),
       overall = F,col.square = "#363886",lwd=2,squaresize = 0.4,xlim=c(0.9,1.6),
       weight.study = "same",col.inside="black",ref = 1,plotwidth="3.5cm")
dev.off()

pdf("multi 14 VS ssgsea CombineSeperate trainC-index.pdf",width=7,height = 7)
Data<-as.data.frame(Data[order(Data[,3],decreasing = T),])
MultiCname<-rownames(Data)
meta.i<-metagen(C.index,se.C.index,sm="C-index",data = Data,
                studlab = rownames(Data))
forest(meta.i, leftcols = c("studlab"),leftlabs = c("Study Element"),
       rightcols = c("effect.ci"),fontsize = 13,just.studlab = "center",
       overall = F,col.square = "#363886",lwd=2,xlim = c(0.4,0.7),
       at=c(0.4,0.5,0.6,0.7),plotwidth="3.5cm",squaresize = 0.4,
       weight.study = "same",col.inside="black",ref = 0.5)
dev.off()
#####10.4 heatmap#####
#####uni#####
HR_C<-unicoxForestData(selfEvalResult)[[2]]
k<-HR_C$lnHR
o<-c()
for(i in 1:length(uni_predictive_ability)){
  o<-rbind(o,uni_predictive_ability[[i]][,1])
}
Data<-cbind(k,o)
rownames(Data)<-rownames(HR_C)
colnames(Data)[c(1,15)]<-c("-FIPI","lu14")
bk<-c(min(Data),seq(-0.5,-0.01,length.out = 49),
      0,seq(0.01,0.5,length.out = 49),max(Data))
color<-c("#7A7ABC",colorRampPalette(c("#7A7ABC","white"))(48),"white",
         colorRampPalette(c("white","firebrick3"))(49),"firebrick3")
Data<-t(Data)
Data<-Data[uniHRname,]
pdf("uni-lnHR-heatmap.pdf",width=6,height=4.5)
pheatmap(Data,cluster_cols = F,cluster_rows = F,
         color=color,breaks=bk,border_color="white",
         cellwidth = 15,cellheight = 15)
dev.off()


k<-HR_C$C_index
o<-c()
for(i in 1:length(uni_predictive_ability)){
  o<-rbind(o,uni_predictive_ability[[i]][,3])
}
Data<-cbind(k,o)
rownames(Data)<-rownames(HR_C)
colnames(Data)[c(1,15)]<-c("ssgsea","lu14")
bk<-c(seq(0.3,0.49,length.out=20),0.5,seq(0.51,0.7,length.out=20),max(Data))
color<-c(colorRampPalette(c("#7A7ABC","white"))(19),"white",
         colorRampPalette(c("white","firebrick3"))(20),"firebrick3")
Data<-t(Data)
Data<-Data[uniCname,]
pdf("uni-C-heatmap.pdf",width=6,height=4.5)
pheatmap(Data,cluster_cols = F,cluster_rows = F,
         color=color,breaks=bk,border_color="white",
         cellwidth = 15,cellheight = 15)
dev.off()


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
pdf("Muti-lnHR-heatmap.pdf",width=6,height=4.5)
pheatmap(Data,cluster_cols = F,cluster_rows = F,
         color=color,breaks=bk,border_color="white",
         cellwidth = 15,cellheight = 15)
dev.off()


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
pdf("Multi-C-heatmap.pdf",width=6,height=5)
pheatmap(Data,cluster_cols = F,cluster_rows = F,
         color=color,breaks=bk,border_color="white",
         cellwidth = 15,cellheight = 15)
dev.off()

#####12.所有数据集整体比较#####
#####12.1ssgsea#####
clinicalMatrix<-getSurvivalMatrix(clinical,"global","Patch et al.")
global<-globalEval(scale(score.all),"global",clinicalMatrix,1,2)
global_Value_se<-data.frame(global[2],global[3],
                            global[9],global[10])

#####12.2 others####
clinicalMatrix<-getSurvivalMatrix(clinical,"global","Patch et al.")
clinicalMatrix[which("Alive"==clinicalMatrix[,2]),2]<-0
clinicalMatrix[which("Dead"==clinicalMatrix[,2]),2]<-1
ForteenPerfo<-c()
input_files<-names(datasets_nonscale)
model_files<-names(model.coefs1)
for(i in 1:length(model.coefs1)){
  modelGlobalScore.i<-c()
  for(j in 1:length(input_files)){
    modelGlobalScore.i<-c(modelGlobalScore.i,
                          modelScore(model_files[i],input_files[j]))
  }
  name<-names(modelGlobalScore.i)
  modelGlobalScore.i<-scale(modelGlobalScore.i)
  if(!is.null(modelGlobalScore.i)){
    globalData<-data.frame(time=as.numeric(clinicalMatrix[match(name,rownames(clinicalMatrix)),1]),
                           status=as.numeric(clinicalMatrix[match(name,rownames(clinicalMatrix)),2]),
                           model_score=modelGlobalScore.i)
    bb<-survConcordance(Surv(time,status)~model_score,globalData)#一致性
    cc<-coxph(Surv(time,status)~model_score,globalData)#HR
    ForteenPerfo<-rbind(ForteenPerfo,c(summary(cc)$coefficients[1,c(1,3)],bb$concordance,bb$std.err))
    rownames(ForteenPerfo)[nrow(ForteenPerfo)]<-names(model.coefs1)[i]
  }
}
colnames(ForteenPerfo)<-c("lnHR","selnHR","C-index","se C-index")
write.csv(ForteenPerfo,"ForteenPerfo.csv")
colnames(ForteenPerfo)<-names(global_Value_se)<-c("lnHR","selnHR","c.index","se.c.index")
Data<-rbind(ForteenPerfo,global_Value_se)
rownames(Data)<-c(rownames(ForteenPerfo),"-FIPI")
pdf("14 VS ssgsea globalTrainC-index.pdf",width=7,height = 7)
Data<-as.data.frame(Data[order(Data[,3],decreasing = T),])
meta.i<-metagen(c.index,se.c.index,sm="C-index",data = Data,
                studlab = rownames(Data))
forest(meta.i, leftcols = c("studlab"),leftlabs = c("Study Element"),
       rightcols = c("effect.ci"),fontsize = 13,just.studlab = "center",
       overall = F,col.square = "#363886",lwd=2,xlim = c(0.4,0.7),
       weight.study = "same",col.inside="black",ref = 0.5)
dev.off()

pdf("14 VS ssgsea globalTrain HR.pdf",width=7,height = 7)
Data<-as.data.frame(Data[order(Data[,1],decreasing = T),])
meta.i<-metagen(lnHR,selnHR,sm="HR",data = Data,
                studlab = rownames(Data))
forest(meta.i, leftcols = c("studlab"),leftlabs = c("Study Element"),
       rightcols = c("effect.ci"),fontsize = 13,just.studlab = "center",
       overall = F,col.square = "#363886",lwd=2,xlim = c(0.9,7),
       weight.study = "same",col.inside="black",ref = 1)
dev.off()
#####13.泡泡图######
setwd("G:\\haodapeng-data\\deleteScatter\\think1")
source('G:/haodapeng-data/deleteScatter/code/Bubble.R')
source('G:/haodapeng-data/deleteScatter/code/Bubble_row_colName.R')
source('G:/haodapeng-data/deleteScatter/code/bubblePlot.R')
library(ggplot2)
#####13.1 相关性泡泡图######
#####13.2 cox gene (term扰动)泡泡图#####
tri_semsim_term<-read.csv("tri_semsim_term.csv",header=T,row.names = 1,as.is=T)
tri_semsim_term_p<-read.csv("tri_semsim_term_p.csv",header=T,row.names = 1,as.is=T)
colnames(tri_semsim_term)<-rownames(tri_semsim_term)
colnames(tri_semsim_term_p)<-rownames(tri_semsim_term_p)
bubblePlot(tri_semsim_term,tri_semsim_term_p,
           "cox sim based on Term permutation",
           "cox sim based on Term permutation.pdf")
#####13.3 cox gene(jiyin扰动)泡泡图#####
load("diff_coxgene.Rdata")
semSimilirity<-read.csv("semSimilirity.csv",header=T,row.names = 1,as.is=T)
semSimilirity_p<-read.csv("semSimilirity_p.csv",header=T,row.names = 1,as.is=T)
colnames(semSimilirity)<-colnames(semSimilirity_p)<-rownames(semSimilirity)<-rownames(semSimilirity_p)<-names(diff.coxgene)
bubblePlot(semSimilirity,semSimilirity_p,
           "cox sim based on gene permutation",
           "cox sim based on gene permutation.pdf")
#####13.4 cox gene(表达相关性)泡泡图#####
load("diff_coxgene.Rdata")
coxGeneSim<-read.csv("coxGeneSim.csv",header=T,row.names = 1,as.is=T)
coxGeneSimP<-read.csv("coxGeneSimP.csv",header=T,row.names = 1,as.is=T)
colnames(coxGeneSim)<-colnames(coxGeneSimP)<-rownames(coxGeneSim)<-rownames(coxGeneSimP)<-names(diff.coxgene)
bubblePlot(coxGeneSim,coxGeneSimP,
           "cox sim based on exp correlation",
           "cox sim based on exp correlation.pdf")
#####13.5 cox gene(jiyin and 表达相关性)泡泡图####
load("diff_coxgene.Rdata")
semSimilirity<-read.csv("semSimilirity.csv",header=T,row.names = 1,as.is=T)
semSimilirity_p<-read.csv("semSimilirity_p.csv",header=T,row.names = 1,as.is=T)
colnames(semSimilirity)<-colnames(semSimilirity_p)<-rownames(semSimilirity)<-rownames(semSimilirity_p)<-names(diff.coxgene)
coxGeneSim<-read.csv("coxGeneSim.csv",header=T,row.names = 1,as.is=T)
coxGeneSimP<-read.csv("coxGeneSimP.csv",header=T,row.names = 1,as.is=T)
colnames(coxGeneSim)<-colnames(coxGeneSimP)<-rownames(coxGeneSim)<-rownames(coxGeneSimP)<-names(diff.coxgene)
# sim<-semSimilirity
# sim[upper.tri(sim)]<-coxGeneSim[upper.tri(coxGeneSim)]
# simP<-semSimilirity_p
# simP[upper.tri(simP)]<-coxGeneSimP[upper.tri(coxGeneSimP)]
# bubblePlot(sim,simP,
#            "cox sim based on exp correlation and gene permutaion",
#            "cox sim based on exp correlation and gene permutaion.pdf")


sim<-semSimilirity
sim[upper.tri(sim)]<-NA
simP<-semSimilirity_p
simP[upper.tri(simP)]<-NA
bubblePlot(sim,simP,
           "plot",
           "plot.pdf")


sim<-coxGeneSim
sim[lower.tri(sim)]<-NA
simP<-coxGeneSimP
simP[lower.tri(simP)]<-NA
bubblePlot(sim,simP,
           "plot1",
           "plot1.pdf")
#####13.6 cox gene (内容)泡泡图#####
jaccardMatrix<-read.csv("jaccardMatrix.csv",header=T,row.names = 1,as.is=T)
jaccardMatrix_p<-read.csv("jaccardMatrix_p.csv",header=T,row.names = 1,as.is=T)
colnames(jaccardMatrix)<-rownames(jaccardMatrix)
colnames(jaccardMatrix_p)<-rownames(jaccardMatrix_p)
bubblePlot(jaccardMatrix,jaccardMatrix_p,
           "cox sim based on jaccard",
           "cox sim based on jaccard.pdf")
#####13.7 cox gene (内容)箱式#####
jaccardMatrix<-read.csv("jaccardMatrix.csv",header=T,row.names = 1,as.is=T)
#jaccardMatrix_p<-read.csv("jaccardMatrix_p.csv",header=T,row.names = 1,as.is=T)
#colnames(jaccardMatrix_p)<-rownames(jaccardMatrix_p)
data<-data.frame(class=rep(rownames(jaccardMatrix),each=15),
                 value=unlist(jaccardMatrix)[!is.na(unlist(jaccardMatrix))])
ggplot(data,aes(class,value))+
  geom_boxplot(fill="#EFC000",outlier.shape = NA)+
  geom_jitter(shape=21,fill="#EFC000",size=1.7,width=0.3,height=0)+
  theme_bw()+
  labs(y="Jaccard Index")+
  theme(panel.border = element_blank(),panel.grid=element_blank(),
      axis.line = element_line(colour = "black"),
      axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1,
                                 color = "black",size=14),
      axis.text.y=element_text(color="black",size=14),
      plot.title = element_text(colour = "black", face = "bold", 
                                size = 14, vjust = 1,hjust = 0.5),
      axis.title.x = element_blank(),
      axis.title.y=element_text(color="black",size=16),
      legend.position ="none")
ggsave("cox sim based on jaccard.pdf",width=6,height=4)

#####14网络图数据整理######
setwd("G:\\haodapeng-data\\deleteScatter\\think1")
library("org.Hs.eg.db")
enrih<-read.csv("enrich_meta.csv",header=T,as.is=T)
enrih1<-matrix(0,nrow = nrow(enrih),ncol=1)
a<-apply(enrih,1,function(x){paste(x[2],x[3],sep = "~")})
enrih1<-cbind(enrih1,enrih[,3],a)
b<-apply(enrih,1,function(x){
  m<-as.numeric(unlist(strsplit(x[4],"/")))
  return(m[1]/m[2])
})
x<-org.Hs.eg.db
c<-apply(enrih,1,function(y){
  i1<-unlist(strsplit(unlist(y[9]),"/",fixed = T))
  entrez<-select(x,keys=i1,columns = "ENTREZID",keytype = "SYMBOL")[[2]]
  return(paste(entrez,"",collapse = ", ",sep = ""))
})
enrih1<-enrih1[,-1]
enrih1<-cbind(enrih1,enrih[,10],b,enrih[,8],c,58,NA,16672,NA,1,1,enrih[,8])
colnames(enrih1)<-c("Category","Term","Count","%","PValue",
                    "Genes","List Total",	"Pop Hits",	"Pop Total",
                    "Fold Enrichment","Bonferroni","Benjamini","FDR")
write.table(enrih1,"enrih1.txt",row.names=F,sep="\t",quote = F)
#####15模型gene注释信息######
#用R 3.6.0版本
library("biomaRt")
gene_enrich_68<-read.csv("gene_enrich_68.csv",header = T,row.name=1,as.is=T)
rownames(gene_enrich_68)[65]<-"SPART"
gene_enrich_68$Class<-ifelse(gene_enrich_68$Class=="Protect","Protective","Risky")
meta.infor<-data.frame(paste(gene_enrich_68$HR,"[",gene_enrich_68$CI,"]",sep=""),
                       gene_enrich_68[,3:6])

ensembl <- useMart("ensembl") 
ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl) 
a<-getBM(attributes = c("hgnc_symbol","ensembl_gene_id","entrezgene","chromosome_name",
                        "start_position",
                     "end_position","strand"),
      filters=c('hgnc_symbol'), 
      values=gene, 
      mart=ensembl)
a<-a[!grepl("C",a$chromosome_name),]
a<-a[match(gene,a$hgnc_symbol),]
a$strand<-ifelse(a$strand==1,"(+)","(-)")
a$chromosome_name<-paste("Chr ",a$chromosome_name,sep="")
a$chroPo<-paste(paste(a$chromosome_name,paste(prettyNum(a$start_position,big.mark = ","),
                                              prettyNum(a$end_position,big.mark = ","),sep="-"),sep=":"),
                a$strand,sep=" ")
b<-data.frame(a$hgnc_symbol,a$ensembl_gene_id,a$entrezgene,a$chroPo)
annotationOFModelGene<-cbind(b,meta.infor)
colnames(annotationOFModelGene)<-c("Gene symbol","Ensembl id","Entrez id",
                                   "Chromosome position","HR(95% CI)","p-value",
                                   "FDR","Class","Weight")
write.csv(annotationOFModelGene,"annotationOFModelGene.csv",row.names = F)

#####16铂治疗患者相关分析#####
#####16.1提取化疗史样本的化疗反应、PFS/RFS#####
setwd("G:\\haodapeng-data\\deleteScatter\\think1")
load("clinical.Rdata")
load("ssgsea_out.Rdata")
load("ssgsea_out_test.Rdata")
load("ssgsea_out_all.Rdata")
load("ssgsea.out.independ.Rdata")
load("datasets_scale.Rdata")
clinical.chemotherapy<-list()
#####train:GSE32062(9th)废弃（无化疗反应数据）,Dressman(13th)####
{#放弃
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
#####test:GSE30161(8th), GSE32063(10th)废弃（无化疗反应数据）####
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
#a$fufa.status<-ifelse(a$fufa.status=="recurrence",1,0)
#po<-unique(c(which(is.na(a$fufa.status)),which(a$chemo.response=="unknown")))
#a<-a[-c(12,15,18),]
clinical.chemotherapy[["GSE30161"]]<-a

{#已放弃
  ctestpo<-10
  preclinical<-t(apply(as.matrix(clinical[[ctestpo]][,31],ncol=1), 1,function(x){unlist(strsplit(x,split = "///"))}))
  rownames(preclinical)<-rownames(clinical[[ctestpo]])
  colnames(preclinical)<-unlist(lapply(strsplit(preclinical[1,],": "), function(x){ifelse(length(x)==2,x[1],x[2])}))
  preclinical<-as.data.frame(t(apply(preclinical,1,function(x){unlist(lapply(strsplit(x,": "), function(x){ifelse(length(x)==2,x[2],x[3])}))})))}
#####test.independ:patch(16th)、TCGA-seq(14th)####
ctest_independpo<-14
a<-data.frame(score=structure(ssgsea.out.independ[[1]][3,],
                              .Names=colnames(ssgsea.out.independ[[1]])),
              fufa.status=clinical[[14]]$RFS,
              fufa.days=clinical[[14]]$RFS.time,
              os.status=clinical[[14]]$vital_status,
              os.days=clinical[[14]]$days_to_death,
              chemo.agent=clinical[[14]]$additional_pharmaceutical_therapy,
              chemo.response=clinical[[14]]$primary_therapy_outcome_success)
#a<-a[-unique(c(which(is.na(a$fufa.status)),which(is.na(a$fufa.days)))),]

#clinical.chemotherapy[["TCGA-seq"]]<-a[25,]
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
save(clinical.chemotherapy,file="clinical_chemotherapy.Rdata")

#####16.2铂治疗各种图#####
library(ggplot2)
source('G:/haodapeng-data/deleteScatter/code/chemoLabel.R')
setwd("G:\\haodapeng-data\\deleteScatter\\think1")
load("clinical_chemotherapy.Rdata")
# criterion.left<-seq(-10500,5000,500)
# criterion.right<-seq(-10000,5500,500)

#####train####
chemo.train<-clinical.chemotherapy[[1]]
chemo.train<-chemo.train[order(chemo.train$score),]
label.train<-c(rep(1:9,each=11),rep(10,9))
# criterion.left<-seq(-6200,1000,length.out=11)[1:10]
# criterion.right<-seq(-6200,1000,length.out=11)[2:11]
# label.train<-apply(chemo.train,1,
#                    function(x){chemo.label(x,criterion.left,criterion.right)})
a<-tapply(chemo.train[,7],factor(label.train,
                                 levels = as.character(sort(unique(as.numeric(label.train))))),function(x){
  return(length(which(x=="completeresponse"))/length(x))
})
Data<-data.frame(score=1:10,
                 CR.propo=as.vector(a))

ggplot(Data,aes(x=score,y=CR.propo))+
  geom_point()+
  #labs(y="Complete response (%)",x=expression(paste("Score"," (",10^3,")")))+
  geom_smooth(method = lm,se=F)

fit<-lm(CR.propo~score,Data)
summary(fit)
#####test I####
criterion.left<-seq(-6200,5100,length.out=11)[1:10]
criterion.right<-seq(-6200,5100,length.out=11)[2:11]
chemo.test1<-clinical.chemotherapy[[2]]
label.test1<-apply(chemo.test1,1,
                   function(x){chemo.label(x,criterion.left,criterion.right)})
a<-tapply(chemo.test1[,7],factor(label.test1,
                                 levels = as.character(sort(unique(as.numeric(label.test1))))),function(x){
                                   return(length(which(x=="CR (completely response)"))/length(x))
                                 })
Data<-data.frame(score=c(1:10)[-9],
                 CR.propo=as.vector(a))

ggplot(Data,aes(x=score,y=CR.propo))+
  geom_point()+
  #labs(y="Complete response (%)",x=expression(paste("Score"," (",10^3,")")))+
  geom_smooth(method = lm,se=F)
fit<-lm(CR.propo~score,Data)
summary(fit)
#####test II####
criterion.left<-seq(-11000,3200,length.out=11)[1:10]
criterion.right<-seq(-11000,3200,length.out=11)[2:11]
chemo.test2<-rbind(clinical.chemotherapy[[3]],clinical.chemotherapy[[4]])
label.test2<-apply(chemo.test2,1,
                   function(x){chemo.label(x,criterion.left,criterion.right)})
a<-tapply(chemo.test2[,7],factor(label.test2,
                                 levels = as.character(sort(unique(as.numeric(label.test2))))),function(x){
                                   return(length(c(which(x=="CR (completely response)"),
                                                   which(x=="Sensitive")))/length(x))
                                 })
Data<-data.frame(score=1:10,
                 CR.propo=as.vector(a))
ggplot(Data,aes(x=score,y=CR.propo))+
  geom_point()+
  #labs(y="Complete response (%)",x=expression(paste("Score"," (",10^3,")")))+
  geom_smooth(method = lm,se=F)
fit<-lm(CR.propo~score,Data)
summary(fit)
#####整合到一起####
write.csv(clinical.chemotherapy[[1]],"clinical.chemotherapy.1.csv")
write.csv(clinical.chemotherapy[[2]],"clinical.chemotherapy.2.csv")
write.csv(clinical.chemotherapy[[3]],"clinical.chemotherapy.3.csv")
write.csv(clinical.chemotherapy[[4]],"clinical.chemotherapy.4.csv")
chemo1<-read.csv("clinical.chemotherapy.1.csv",header=T,row.names = 1,as.is=T)
chemo2<-read.csv("clinical.chemotherapy.2.csv",header=T,row.names = 1,as.is=T)
chemo3<-read.csv("clinical.chemotherapy.3.csv",header=T,row.names = 1,as.is=T)
chemo4<-read.csv("clinical.chemotherapy.4.csv",header=T,row.names = 1,as.is=T)
chemo<-rbind(rbind(rbind(chemo1,chemo2),chemo3),chemo4)
#####按分数等分####
criterion.left<-seq(-11000,5100,length.out=11)[1:10]
criterion.right<-seq(-11000,5100,length.out=11)[2:11]
label<-apply(chemo,1,function(x){chemo.label(x,criterion.left,criterion.right)})
a<-tapply(chemo[,7],factor(label,
                           levels = as.character(sort(unique(as.numeric(label.test2))))),
          function(x){return(length(c(which(x=="CR (completely response)"),
                                      which(x=="Sensitive"),
                                      which(x=="Complete Remission/Response"),
                                      which(x=="completeresponse")))/length(x))
          })
#####按样本数等分####
#散点图-global
chemo<-chemo[order(chemo$score),]
write.csv(chemo,"chemo.csv")
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
#散点图-train
chemo1<-chemo1[order(chemo1$score),]
label<-c(rep(1:9,each=11),rep(10,9))
a<-c()
for(i in 1:10){
  x<-c()
  if(i!=10){
    x<-chemo1[(11*(i-1)+1):(11*i),7]
  }else{
    x<-chemo1[100:108,7]
  }
  a<-c(a,length(c(which(x=="CR (completely response)"),
                  which(x=="Sensitive"),
                  which(x=="Complete Remission/Response"),
                  which(x=="completeresponse")))/length(x))
}
Data<-data.frame(phase=1:10,
                 CR.propo=a)

ggplot(Data,aes(x=phase,y=CR.propo))+
  geom_point()+
  #labs(y="Complete response (%)",x=expression(paste("Score"," (",10^3,")")))+
  geom_smooth(method = lm,se=F)

fit<-lm(CR.propo~phase,Data)
summary(fit)

#散点图-test1
chemo2<-chemo2[order(chemo2$score),]
label<-c(rep(1:9,each=4),rep(10,6))
a<-c()
for(i in 1:10){
  x<-c()
  if(i!=10){
    x<-chemo2[(4*(i-1)+1):(4*i),7]
  }else{
    x<-chemo2[37:442,7]
  }
  a<-c(a,length(c(which(x=="CR (completely response)"),
                  which(x=="Sensitive"),
                  which(x=="Complete Remission/Response"),
                  which(x=="completeresponse")))/length(x))
}
Data<-data.frame(phase=1:10,
                 CR.propo=a)

ggplot(Data,aes(x=phase,y=CR.propo))+
  geom_point()+
  #labs(y="Complete response (%)",x=expression(paste("Score"," (",10^3,")")))+
  geom_smooth(method = lm,se=F)

fit<-lm(CR.propo~phase,Data)
summary(fit)
#散点图-test2
chemo34<-rbind(chemo3,chemo4)
chemo34<-chemo34[order(chemo34$score),]
label<-rep(1:10,each=28)
a<-c()
for(i in 1:10){
  x<-chemo34[(28*(i-1)+1):(28*i),7]
  a<-c(a,length(c(which(x=="CR (completely response)"),
                  which(x=="Sensitive"),
                  which(x=="Complete Remission/Response"),
                  which(x=="completeresponse")))/length(x))
}
Data<-data.frame(phase=1:10,
                 CR.propo=a)

ggplot(Data,aes(x=phase,y=CR.propo))+
  geom_point()+
  #labs(y="Complete response (%)",x=expression(paste("Score"," (",10^3,")")))+
  geom_smooth(method = lm,se=F)

fit<-lm(CR.propo~phase,Data)
summary(fit)

#散点图-test总
chemo234<-rbind(chemo2,chemo3,chemo4)
chemo234<-chemo234[order(chemo234$score),]
label<-c(rep(1:9,each=32),rep(10,34))
a<-c()
for(i in 1:10){
  x<-c()
  if(i!=10){
    x<-chemo234[(32*(i-1)+1):(32*i),7]
  }else{
    x<-chemo234[289:322,7]
  }
  a<-c(a,length(c(which(x=="CR (completely response)"),
                  which(x=="Sensitive"),
                  which(x=="Complete Remission/Response"),
                  which(x=="completeresponse")))/length(x))
}
Data<-data.frame(phase=1:10,
                 CR.propo=a)

ggplot(Data,aes(x=phase,y=CR.propo))+
  geom_point()+
  #labs(y="Complete response (%)",x=expression(paste("Score"," (",10^3,")")))+
  geom_smooth(method = lm,se=F)

fit<-lm(CR.propo~phase,Data)
summary(fit)
###ROC####
#global
{
  library(pROC)
  chemo$chemo.label<-ifelse(chemo$chemo.response=="CR (completely response)"|
                              chemo$chemo.response=="Sensitive"|
                              chemo$chemo.response=="Complete Remission/Response"|
                              chemo$chemo.response=="completeresponse"|
                              chemo$chemo.response=="Partial Remission/Response"|
                              chemo$chemo.response=="PR (partial response)",1,0)
  roc.result<-roc(chemo.label~score,chemo,direction=">")
  plot(roc.result)
}

###train,test1,test2
{
  chemo1$chemo.label<-ifelse(chemo1$chemo.response=="CR (completely response)"|
                               chemo1$chemo.response=="Sensitive"|
                               chemo1$chemo.response=="Complete Remission/Response"|
                               chemo1$chemo.response=="completeresponse"|
                               chemo1$chemo.response=="Partial Remission/Response"|
                               chemo1$chemo.response=="PR (partial response)",1,0)
  roc(chemo.label~score,chemo1,direction=">")
  
  chemo2$chemo.label<-ifelse(chemo2$chemo.response=="CR (completely response)"|
                               chemo2$chemo.response=="Sensitive"|
                               chemo2$chemo.response=="Complete Remission/Response"|
                               chemo2$chemo.response=="completeresponse"|
                               chemo2$chemo.response=="Partial Remission/Response"|
                               chemo2$chemo.response=="PR (partial response)",1,0)
  roc(chemo.label~score,chemo2,direction=">")
  
  chemo34<-rbind(chemo3,chemo4)
  chemo34$chemo.label<-ifelse(chemo34$chemo.response=="CR (completely response)"|
                                chemo34$chemo.response=="Sensitive"|
                                chemo34$chemo.response=="Complete Remission/Response"|
                                chemo34$chemo.response=="completeresponse"|
                                chemo34$chemo.response=="Partial Remission/Response"|
                                chemo34$chemo.response=="PR (partial response)",1,0)
  roc(chemo.label~score,chemo34,direction=">")
}

###array vs seq:CR+PR
{
  chemo.array<-rbind(rbind(chemo1,chemo2),chemo4)
  chemo.array$chemo.label<-ifelse(chemo.array$chemo.response=="CR (completely response)"|
                                    chemo.array$chemo.response=="Sensitive"|
                                    chemo.array$chemo.response=="Complete Remission/Response"|
                                    chemo.array$chemo.response=="completeresponse"|
                                    chemo.array$chemo.response=="Partial Remission/Response"|
                                    chemo.array$chemo.response=="PR (partial response)",1,0)
  roc(chemo.label~score,chemo.array,direction=">")
  
  chemo.seq<-chemo3
  chemo.seq$chemo.label<-ifelse(chemo.seq$chemo.response=="CR (completely response)"|
                                  chemo.seq$chemo.response=="Sensitive"|
                                  chemo.seq$chemo.response=="Complete Remission/Response"|
                                  chemo.seq$chemo.response=="completeresponse"|
                                  chemo.seq$chemo.response=="Partial Remission/Response"|
                                  chemo.seq$chemo.response=="PR (partial response)",1,0)
  roc(chemo.label~score,chemo.seq,direction=">")
}
###array vs seq:CR
{
  chemo.array.1<-rbind(rbind(chemo1,chemo2),chemo4)
  chemo.array.1$chemo.label<-ifelse(chemo.array.1$chemo.response=="CR (completely response)"|
                                      chemo.array.1$chemo.response=="Sensitive"|
                                      chemo.array.1$chemo.response=="Complete Remission/Response"|
                                      chemo.array.1$chemo.response=="completeresponse",1,0)
  roc(chemo.label~score,chemo.array.1,direction=">")
  
  chemo.seq.1<-chemo3
  chemo.seq.1$chemo.label<-ifelse(chemo.seq.1$chemo.response=="CR (completely response)"|
                                    chemo.seq.1$chemo.response=="Sensitive"|
                                    chemo.seq.1$chemo.response=="Complete Remission/Response"|
                                    chemo.seq.1$chemo.response=="completeresponse",1,0)
  roc(chemo.label~score,chemo.seq.1,direction=">")
}
###shengcun#####
#medianscore<-2687.859
library(survival)
library(survminer)
chemo.cr<-chemo[chemo$chemo.label==1,]
medianscore<-median(chemo.cr$score)
chemo.cr$os.status.1<-ifelse(chemo.cr$os.status=="0"|chemo.cr$os.status=="deceased",0,1)
chemo.cr$group<-ifelse(chemo.cr$score>medianscore,"High Risk","Low Risk")
fit<-survfit(Surv(os.days,os.status.1)~group,data=chemo.cr)
diff<-survdiff(Surv(os.days,os.status.1)~group,data=chemo.cr)
ggsurvplot(fit, data = chemo.cr, title="",
           risk.table = TRUE,risk.table.col = "strata",pval = 0.11,
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
###train,test1,test2
{
  chemo1.cr<-chemo1[chemo1$chemo.label==1,]
  medianscore<-median(chemo1.cr$score)
  chemo1.cr$os.status.1<-ifelse(chemo1.cr$os.status=="0"|chemo1.cr$os.status=="deceased",0,1)
  chemo1.cr$group<-ifelse(chemo1.cr$score>medianscore,"High Risk","Low Risk")
  fit<-survfit(Surv(os.days,os.status.1)~group,data=chemo1.cr)
  diff<-survdiff(Surv(os.days,os.status.1)~group,data=chemo1.cr)
  ggsurvplot(fit, data = chemo1.cr, title="",
             risk.table = TRUE,risk.table.col = "strata",pval = 0.11,
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
  
  
  chemo2.cr<-chemo2[chemo2$chemo.label==1,]
  medianscore<-median(chemo2.cr$score)
  chemo2.cr$os.status.1<-ifelse(chemo2.cr$os.status=="0"|chemo2.cr$os.status=="deceased",0,1)
  chemo2.cr$group<-ifelse(chemo2.cr$score>medianscore,"High Risk","Low Risk")
  fit<-survfit(Surv(os.days,os.status.1)~group,data=chemo2.cr)
  diff<-survdiff(Surv(os.days,os.status.1)~group,data=chemo2.cr)
  ggsurvplot(fit, data = chemo2.cr, title="",
             risk.table = TRUE,risk.table.col = "strata",pval = 0.11,
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
  
  
  chemo34.cr<-chemo34[chemo34$chemo.label==1,]
  medianscore<-median(chemo34.cr$score)
  chemo34.cr$os.status.1<-ifelse(chemo34.cr$os.status=="0"|chemo34.cr$os.status=="deceased",0,1)
  chemo34.cr$group<-ifelse(chemo34.cr$score>medianscore,"High Risk","Low Risk")
  fit<-survfit(Surv(os.days,os.status.1)~group,data=chemo34.cr)
  diff<-survdiff(Surv(os.days,os.status.1)~group,data=chemo34.cr)
  ggsurvplot(fit, data = chemo34.cr, title="",
             risk.table = TRUE,risk.table.col = "strata",pval = 0.11,
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
  
}
###array vs seq:CR+PR
{
  chemo.array.cr.pr<-chemo.array[chemo.array$chemo.label==1,]
  medianscore<-median(chemo.array.cr.pr$score)
  chemo.array.cr.pr$os.status.1<-ifelse(chemo.array.cr.pr$os.status=="0"|chemo.array.cr.pr$os.status=="deceased",0,1)
  chemo.array.cr.pr$group<-ifelse(chemo.array.cr.pr$score>medianscore,"High Risk","Low Risk")
  fit<-survfit(Surv(os.days,os.status.1)~group,data=chemo.array.cr.pr)
  diff<-survdiff(Surv(os.days,os.status.1)~group,data=chemo.array.cr.pr)
  ggsurvplot(fit, data = chemo.array.cr.pr, title="array",
             risk.table = TRUE,risk.table.col = "strata",pval = 0.733,
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
  
  chemo.seq.cr.pr<-chemo.seq[chemo.seq$chemo.label==1,]
  medianscore<-median(chemo.seq.cr.pr$score)
  chemo.seq.cr.pr$os.status.1<-ifelse(chemo.seq.cr.pr$os.status=="0"|chemo.seq.cr.pr$os.status=="deceased",0,1)
  chemo.seq.cr.pr$group<-ifelse(chemo.seq.cr.pr$score>medianscore,"High Risk","Low Risk")
  fit<-survfit(Surv(os.days,os.status.1)~group,data=chemo.seq.cr.pr)
  diff<-survdiff(Surv(os.days,os.status.1)~group,data=chemo.seq.cr.pr)
  ggsurvplot(fit, data = chemo.seq.cr.pr, title="seq",
             risk.table = TRUE,risk.table.col = "strata",pval = 0.0493,
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
}
###array vs seq:CR
{
  chemo.array.cr<-chemo.array.1[chemo.array.1$chemo.label==1,]
  medianscore<-median(chemo.array.cr$score)
  chemo.array.cr$os.status.1<-ifelse(chemo.array.cr$os.status=="0"|chemo.array.cr$os.status=="deceased",0,1)
  chemo.array.cr$group<-ifelse(chemo.array.cr$score>medianscore,"High Risk","Low Risk")
  fit<-survfit(Surv(os.days,os.status.1)~group,data=chemo.array.cr)
  diff<-survdiff(Surv(os.days,os.status.1)~group,data=chemo.array.cr)
  ggsurvplot(fit, data = chemo.array.cr, title="array",
             risk.table = TRUE,risk.table.col = "strata",pval = 0.733,
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
  
  chemo.seq.cr<-chemo.seq.1[chemo.seq.1$chemo.label==1,]
  medianscore<-median(chemo.seq.cr$score)
  chemo.seq.cr$os.status.1<-ifelse(chemo.seq.cr$os.status=="0"|chemo.seq.cr$os.status=="deceased",0,1)
  chemo.seq.cr$group<-ifelse(chemo.seq.cr$score>medianscore,"High Risk","Low Risk")
  fit<-survfit(Surv(os.days,os.status.1)~group,data=chemo.seq.cr)
  diff<-survdiff(Surv(os.days,os.status.1)~group,data=chemo.seq.cr)
  ggsurvplot(fit, data = chemo.seq.cr, title="seq",
             risk.table = TRUE,risk.table.col = "strata",pval = 0.367,
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
}
#####17PFS/RFS区分样本####
library(survival)
library(survminer)
setwd("G:\\haodapeng-data\\deleteScatter\\think1")
load("clinical.Rdata")
load("ssgsea_out.Rdata")
load("ssgsea_out_test.Rdata")
load("ssgsea_out_all.Rdata")
load("ssgsea.out.independ.Rdata")
medianScore<-2687.859
#####17.1所有数据集生存数据类型####
i=14
clinical[[i]][1,]
clinical[[i]]$RFS.time
clinical[[i]]$RFS
table(clinical[[i]]$recurrence_status)
length(clinical[[i]]$recurrence_status)

length(clinical[[i]]$days_to_death)
clinical[[i]]$days_to_death
table(clinical[[i]]$vital_status)


i=10

clinical[[i]]$days_to_tumor_recurrence
clinical[[i]]$recurrence_status

#####17.2 Training 、test、independent test#####
#train

trainData<-data.frame(rbind(clinical[[11]][,c(16,17)],
                            clinical[[12]][,c(16,17)]),
                      score=c(ssgsea.out[[5]][3,],
                            ssgsea.out[[6]][3,]))
#names(clinical)[c(11,12)]
#names(ssgsea.out)[c(5,6)]
trainData<-trainData[-which(is.na(trainData[,1])),]
trainData$group<-ifelse(trainData$score>medianScore,"High Risk","Low Risk")
trainData$recurrence_status<-ifelse(trainData$recurrence_status=="recurrence",1,0)
fit<-survfit(Surv(days_to_tumor_recurrence,recurrence_status)~group,
             data=trainData)
diff<-survdiff(Surv(days_to_tumor_recurrence,recurrence_status)~group,
               data=trainData)
ggsurvplot(fit, data = trainData, title="Training dataset(n = 377)",
           risk.table = TRUE,risk.table.col = "strata",
           pval =  format(1-pchisq(diff$chisq,1),digits = 3),
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
#test
test1Data<-data.frame(rbind(clinical[[4]][,c(16,17)],
                            clinical[[6]][,c(16,17)],
                            clinical[[8]][,c(16,17)]),
                      score=c(ssgsea.out.test[[2]][3,],
                              ssgsea.out.test[[4]][3,],
                              ssgsea.out.test[[5]][3,]))
#names(clinical)[c(4,6,8)]
#names(ssgsea.out.test)[c(2,4,5)]
test1Data<-test1Data[-which(is.na(test1Data[,2])),]
test1Data$group<-ifelse(test1Data$score>medianScore,"High Risk","Low Risk")
test1Data$recurrence_status<-ifelse(test1Data$recurrence_status=="recurrence",1,0)
fit<-survfit(Surv(days_to_tumor_recurrence,recurrence_status)~group,
             data=test1Data)
diff<-survdiff(Surv(days_to_tumor_recurrence,recurrence_status)~group,
               data=test1Data)
ggsurvplot(fit, data = test1Data, title="Test dataset(n = 186)",
           risk.table = TRUE,risk.table.col = "strata",
           pval =  format(1-pchisq(diff$chisq,1),digits = 3),
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
#independent test
clinical[[16]][,6]<-ifelse(clinical[[16]][,6]=="relapse",1,0)
test2Data<-data.frame(days_to_tumor_recurrence=c(clinical[[14]][,16],clinical[[16]][,11]),
                      recurrence_status=c(clinical[[14]][,17],clinical[[16]][,6]),
                      score=c(ssgsea.out.independ[[1]][3,],
                              ssgsea.out.independ[[2]][3,]))
#names(clinical)[c(14,16)]
#names(ssgsea.out.independ)[c(1,2)]
test2Data<-test2Data[-which(is.na(test2Data[,1])),]
test2Data$group<-ifelse(test2Data$score>medianScore,"High Risk","Low Risk")

fit<-survfit(Surv(days_to_tumor_recurrence,recurrence_status)~group,
             data=test2Data)
diff<-survdiff(Surv(days_to_tumor_recurrence,recurrence_status)~group,
               data=test2Data)
ggsurvplot(fit, data = test2Data, title="Test dataset(n = 186)",
           risk.table = TRUE,risk.table.col = "strata",
           pval =  format(1-pchisq(diff$chisq,1),digits = 3),
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
#####17.3 array、seq#####
#array
arrayData<-data.frame(rbind(clinical[[11]][,c(16,17)],
                            clinical[[12]][,c(16,17)],
                            clinical[[4]][,c(16,17)],
                            clinical[[6]][,c(16,17)],
                            clinical[[8]][,c(16,17)]),
                      score=c(ssgsea.out[[5]][3,],
                              ssgsea.out[[6]][3,],
                              ssgsea.out.test[[2]][3,],
                              ssgsea.out.test[[4]][3,],
                              ssgsea.out.test[[5]][3,]))
#names(clinical)[c(11,12)]
#names(ssgsea.out)[c(5,6)]
arrayData<-arrayData[-which(is.na(arrayData[,1])),]
arrayData<-arrayData[-which(is.na(arrayData[,2])),]
arrayData$group<-ifelse(arrayData$score>medianScore,"High Risk","Low Risk")
arrayData$recurrence_status<-ifelse(arrayData$recurrence_status=="recurrence",1,0)
fit<-survfit(Surv(days_to_tumor_recurrence,recurrence_status)~group,
             data=arrayData)
diff<-survdiff(Surv(days_to_tumor_recurrence,recurrence_status)~group,
               data=arrayData)
ggsurvplot(fit, data = arrayData, title="Array dataset(n = 563)",
           risk.table = TRUE,risk.table.col = "strata",
           pval =  format(1-pchisq(diff$chisq,1),digits = 3),
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
#####17.4 所有数据集合并到一起#####
globalData<-rbind(trainData,test1Data,test2Data)
fit<-survfit(Surv(days_to_tumor_recurrence,recurrence_status)~group,
             data=globalData)
diff<-survdiff(Surv(days_to_tumor_recurrence,recurrence_status)~group,
               data=globalData)
ggsurvplot(fit, data = globalData, title="Array dataset(n = 667)",
           risk.table = TRUE,risk.table.col = "strata",
           pval =  format(1-pchisq(diff$chisq,1),digits = 3),
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

#####18同行评议#####
#####18.1常规生存分析signature构建#####
source('G:/haodapeng-data/deleteScatter/code/getSurvivalMatrix.R')
source('G:/haodapeng-data/deleteScatter/code/getSurvivalMatrix1.R')
source('G:/haodapeng-data/deleteScatter/code/mergeexp.R')
library(survival)
library(survminer)
setwd("G:/haodapeng-data/deleteScatter/think1")
load("datasets_scale.Rdata")
load("datasets_nonscale.Rdata")
load("clinical.Rdata")
trainPo<-which(sapply(datasets_scale,ncol)>100)[-8]
testPo<-which(sapply(datasets_scale,ncol)<100)[-8]
independPo<-c(14,16)
oriexp<-datasets_nonscale
#oriexp<-datasets_scale
#####18.1.1training#####
#生存数据
clinicalTrain<-clinical[trainPo]
clinicalMatrix<-getSurvivalMatrix(clinicalTrain,"global","Patch et al.")
clinicalMatrix[,2]<-ifelse(clinicalMatrix[,2]=="Dead",1,0)
# #clinicalMatrix<-as.data.frame(t(apply(clinicalMatrix,1,as.numeric)))
# colnames(clinicalMatrix)<-c("time","status")
#表达谱数据
#uniongene<-unique(unlist(sapply(list,rownames)))
list<-oriexp[trainPo]
lll<-table(unlist(sapply(oriexp,rownames)))
uniongene<-names(lll)[lll==16]
trainexp<-matrix(NA,nrow=3615,ncol=1190)
rownames(trainexp)<-uniongene
colnames(trainexp)<-unlist(sapply(list, colnames))
for(i in 1:length(list)){
  rowpo.list<-match(rownames(trainexp),rownames(list[[i]]))
  colpo<-match(colnames(list[[i]]),colnames(trainexp))
  trainexp[,colpo]<-list[[i]][rowpo.list,]
}
#phe<-as.data.frame(clinicalMatrix)
a<-c()
for(i in 1:nrow(trainexp)){
  # po<-which(!is.na(trainexp[i,]))
  # exp.i<-trainexp[i,po]
  # phe.i<-as.data.frame(clinicalMatrix[po,])
  exp.i<-trainexp[i,]
  #gene.i <- coxph(Surv(as.numeric(phe$time),as.numeric(phe$status))~ exp.i)
  gene.i <- coxph(Surv(as.numeric(clinicalMatrix[,1]),
                       as.numeric(clinicalMatrix[,2]))~ exp.i)
  b<-summary(gene.i)
  a<-rbind(a,c(b$coefficients,b$conf.int[1,2:4]))
}
colnames(a)<-c("coef","exp(coef)","se(coef)","z",
               "Pr(>|z|)","exp(-coef)",
               "lower .95","upper .95")
rownames(a)<-uniongene
write.csv(a,"traincoxgene.csv",quote=F)
#a<-read.csv("traincoxgene.csv",as.is=T,row.names=1,header=T)
a<-cbind(a,p.adjust(a[,5],method="fdr"))
candiGene<-rownames(a)[a[,9]<0.01]
# candiGene<-rownames(a)[a[,9]<0.05]
trainexp1<-as.data.frame(t(trainexp[match(candiGene,rownames(trainexp)),]))

multicox<-coxph(Surv(as.numeric(clinicalMatrix[,1]),
                     as.numeric(clinicalMatrix[,2]))~.,data = trainexp1)
multicox.sum<-summary(multicox)
# length(which(multicox.sum$coefficients[,5]<0.05))
finalmodelgene<-rownames(multicox.sum$coefficients)[which(multicox.sum$coefficients[,5]<0.05)]
trainexp2<-as.data.frame(trainexp1[,match(finalmodelgene,colnames(trainexp1))])
multicox1<-coxph(Surv(as.numeric(clinicalMatrix[,1]),
                      as.numeric(clinicalMatrix[,2]))~.,data = trainexp2)
multicox1.sum<-summary(multicox1)
coxmodel<-structure(multicox1.sum$coefficients[,1],
                    names=rownames(multicox1.sum$coefficients))
save(coxmodel,file="coxmodel.Rdata")

#score
Data<-trainexp[match(names(coxmodel),rownames(trainexp)),]
score<-structure(apply(Data,2,function(x){return(sum((x*coxmodel)))}),
                 names=colnames(Data))
MEDIAN<-median(score)#-0.002533384  
group<-ifelse(score>=MEDIAN,"HighRisk","LowRisk")
survivalMatrix<-as.data.frame(cbind(clinicalMatrix,score,group))
fit<-survfit(Surv(as.numeric(as.character(survivalMatrix$time)),
                  as.numeric(as.character(survivalMatrix$status)))~group,
             data=survivalMatrix)
diff<-survdiff(Surv(as.numeric(as.character(survivalMatrix$time)),
                    as.numeric(as.character(survivalMatrix$status)))~group,
               data=survivalMatrix)
ggsurvplot(fit, data = survivalMatrix, title="train(n = 1190)",
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
#####18.1.2test1#####
clinicalTest<-clinical[testPo]
clinicalMatrix<-getSurvivalMatrix(clinicalTest,"global","Patch et al.")
clinicalMatrix[,2]<-ifelse(clinicalMatrix[,2]=="Dead",1,0)
# clinicalMatrix<-as.data.frame(t(apply(clinicalMatrix,1,as.numeric)))
# colnames(clinicalMatrix)<-c("time","status")
#load("coxmodel.Rdata")
#score
score.test<-c()
for(i in 1:length(oriexp[testPo])){
  x<-oriexp[testPo][[i]]
  if(!any(is.na(match(names(coxmodel),rownames(x))))){
    data.i<-x[match(names(coxmodel),rownames(x)),]
    score.test<-c(score.test,apply(data.i,2,function(x){return(sum((x*coxmodel)))}))
  }else{
    o<-intersect(names(coxmodel),rownames(x))
    data.i<-x[match(o,rownames(x)),]
    coxmodel.i<-coxmodel[match(o,names(coxmodel))]
    score.test<-c(score.test,apply(data.i,2,function(x){return(sum((x*coxmodel.i)))}))
  }
}
#MEDIAN<--0.002533384  
group.test<-ifelse(score.test>=MEDIAN,"HighRisk","LowRisk")
survivalMatrix<-as.data.frame(cbind(clinicalMatrix,score.test,group.test))
fit<-survfit(Surv(as.numeric(as.character(survivalMatrix$time)),
                  as.numeric(as.character(survivalMatrix$status)))~group.test,
             data=survivalMatrix)
diff<-survdiff(Surv(as.numeric(as.character(survivalMatrix$time)),
                    as.numeric(as.character(survivalMatrix$status)))~group.test,
               data=survivalMatrix)
ggsurvplot(fit, data = survivalMatrix, title="test1(n = 438)",
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
#names(score.test)<-unlist(lapply(datasets_scale[testPo], colnames))
#####18.1.3test2#####
clinicalIndepend<-clinical[independPo]#TCGA、patch
clinicalMatrix11<-getSurvivalMatrix1(clinicalIndepend[1],"global","Patch et al.")
clinicalMatrix22<-getSurvivalMatrix(clinicalIndepend[2],"global","Patch et al.")
clinicalMatrix<-rbind(clinicalMatrix11,clinicalMatrix22)
clinicalMatrix[,2]<-ifelse(clinicalMatrix[,2]=="Dead",1,0)
# clinicalMatrix<-as.data.frame(t(apply(clinicalMatrix,1,as.numeric)))
# colnames(clinicalMatrix)<-c("time","status")

score.test1<-c()
for(i in 1:length(oriexp[independPo])){
  x<-oriexp[independPo][[i]]
  if(!any(is.na(match(names(coxmodel),rownames(x))))){
    data.i<-x[match(names(coxmodel),rownames(x)),]
    score.test1<-c(score.test1,apply(data.i,2,function(x){return(sum((x*coxmodel)))}))
  }else{
    o<-intersect(names(coxmodel),rownames(x))
    data.i<-x[match(o,rownames(x)),]
    coxmodel.i<-coxmodel[match(o,names(coxmodel))]
    score.test1<-c(score.test1,apply(data.i,2,function(x){return(sum((x*coxmodel.i)))}))
  }
}
#MEDIAN<--0.002533384
group.test1<-ifelse(score.test1>=MEDIAN,"HighRisk","LowRisk")
survivalMatrix<-as.data.frame(cbind(clinicalMatrix,score.test1,group.test1))
fit<-survfit(Surv(as.numeric(as.character(survivalMatrix$time)),
                  as.numeric(as.character(survivalMatrix$status)))~group.test1,
             data=survivalMatrix)
diff<-survdiff(Surv(as.numeric(as.character(survivalMatrix$time)),
                    as.numeric(as.character(survivalMatrix$status)))~group.test1,
               data=survivalMatrix)
ggsurvplot(fit, data = survivalMatrix, title="test2(n = 354)",
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

#####18.2ssgsea模型与所有模型ROC的比较#####
#####前期准备#####
source('G:/haodapeng-data/deleteScatter/code/signatureScore.R')
source('G:/haodapeng-data/deleteScatter/code/getTotalGene.R')
source('G:/haodapeng-data/deleteScatter/code/ROC_CI_global.R')
source('G:/haodapeng-data/deleteScatter/code/globalEval.R')
source('G:/haodapeng-data/deleteScatter/code/HGClinicalStatus.R')
source('G:/haodapeng-data/deleteScatter/code/plotSurvivalRoc.R')
source('G:/haodapeng-data/deleteScatter/code/sstselfEvalUMContinue.R')
source('G:/haodapeng-data/deleteScatter/code/selfEval.R')
source('G:/haodapeng-data/deleteScatter/code/sstselfEvalUM.R')
source('G:/haodapeng-data/deleteScatter/code/unicoxForestData.R')
source('G:/haodapeng-data/deleteScatter/code/MulticoxForestData.R')
source('G:/haodapeng-data/deleteScatter/code/sstselfEval.R')
source('G:/haodapeng-data/deleteScatter/code/groupQ.R')
source('G:/haodapeng-data/deleteScatter/code/createTable.R')
source('G:/haodapeng-data/deleteScatter/code/ROC_CI.R')
source('G:/haodapeng-data/deleteScatter/think1/combineCoxResult.R')
source('G:/haodapeng-data/deleteScatter/code/getSurvivalMatrix.R')
source('G:/haodapeng-data/deleteScatter/code/plotSurvivalCurve.R')
source('G:/haodapeng-data/deleteScatter/code/getSurvivalMatrix1.R')
source('G:/haodapeng-data/deleteScatter/code/group.R')
source('G:/haodapeng-data/deleteScatter/code/Digits4.R')
source('G:/haodapeng-data/deleteScatter/code/getGlobalExp.R')
source('G:/haodapeng-data/deleteScatter/code/getGlobalScore.R')
source('G:/haodapeng-data/deleteScatter/code/getHeatmapExp.R')
source('G:/haodapeng-data/deleteScatter/code/ROC3_5.R')
source('G:/haodapeng-data/deleteScatter/code/ROC3_5_CI.R')

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
# source('G:/haodapeng-data/deleteScatter/code/timeROC.R')
# source('G:/haodapeng-data/deleteScatter/code/Compute_iid_KM.R')
# source('G:/haodapeng-data/deleteScatter/code/compute_iid_decomposition.R')
# source('G:/haodapeng-data/deleteScatter/code/compute_iid_decomposition_survival.R')
setwd("G:\\haodapeng-data\\deleteScatter\\think1")
load("model_method.RData")
load("model_coefs1.Rdata")
load("datasets_nonscale.Rdata")
load("datasets_scale.Rdata")
load("ssgsea_out.Rdata")
load("ssgsea_out_test.Rdata")
load("ssgsea_out_all.Rdata")
load("ssgsea.out.independ.Rdata")
# load("ssgsea.out.tcga.Rdata")
# load("ssgsea.out.tcga.cel.Rdata")
# load("ssgsea.out.patch.Rdata")
load("clinical.Rdata")
protectModelGene<-as.character(unlist(read.table("protectModelGene.txt",as.is=T)))
riskModelGene<-as.character(unlist(read.table("riskModelGene.txt",as.is=T)))
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
# score.test.tcga<-structure(as.numeric(unlist(sapply(ssgsea.out.tcga,function(x){x[3,]}))),
#                            .Names=as.character(unlist(sapply(ssgsea.out.tcga,colnames))))
score.test<-structure(as.numeric(unlist(sapply(ssgsea.out.test,function(x){x[3,]}))),
                      .Names=as.character(unlist(sapply(ssgsea.out.test,colnames))))
score.all<-structure(as.numeric(unlist(sapply(ssgsea.out.all,function(x){x[3,]}))),
                     .Names=as.character(unlist(sapply(ssgsea.out.all,colnames))))
score.test.independ<-structure(as.numeric(unlist(sapply(ssgsea.out.independ,function(x){x[3,]}))),
                               .Names=as.character(unlist(sapply(ssgsea.out.independ,colnames))))

# pdf("global Trainingsst.pdf",width=6,height=6)
# globalTrain<-globalEval(scale(score),"global Trainingsst",clinicalMatrixTrain,1,2)
# dev.off()


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

model.coef1.cox<-model.coefs1[c(1,2,3,5,6,11)]#6个cox模型
model_files<-names(model.coef1.cox)
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

#ROC3_5(scale(score),clinicalMatrixTrain,1,2)
ppp<-c()
for(i in 1:length(model7score.train)){
  #ppp<-rbind(ppp,ROC3_5(scale(model7score.train[[i]]),clinicalMatrixTrain,1,2))
  ppp<-rbind(ppp,ROC3_5_CI(scale(model7score.train[[i]]),clinicalMatrixTrain,1,2))
  
}
rownames(ppp)<-names(model7score.train)
colnames(ppp)<-c("3year","5year")
allROC[["ssgsea_cox"]][["train"]]<-ppp


data<-data.frame(Model=rep(rownames(ppp),2),
                 AUC=c(ppp[,1],ppp[,2]),
                 Group=c(rep("3 years",7),
                         rep("5 years",7)))
#####test1#####
datasets_scale_test<-datasets_scale[testPo]
datasets_nonscale_test<-datasets_nonscale[testPo]
input_files<-names(datasets_nonscale)[testPo]

model7score.test<-list()
for(i in 1:length(model.coef1)){
  modelGlobalScore.i<-c()
  for(j in 1:length(input_files)){
    modelGlobalScore.i<-c(modelGlobalScore.i,
                          modelScore(model_files[i],input_files[j]))
  }
  model7score.test[[names(model.coef1)[i]]]<-modelGlobalScore.i
  
}
model7score.test[["ssgsea"]]<-score.test

#ROC3_5(scale(score),clinicalMatrixtest,1,2)
ppp<-c()
for(i in 1:length(model7score.test)){
  score.i<-model7score.test[[i]]
  phe.i<-clinicalMatrixTest[match(names(score.i),rownames(clinicalMatrixTest)),]
  ppp<-rbind(ppp,ROC3_5_CI(scale(score.i),phe.i,1,2))
}
rownames(ppp)<-names(model7score.test)
colnames(ppp)<-c("3year","5year")
allROC[["ssgsea_cox"]][["test"]]<-ppp
#####test2#####
datasets_scale_independ<-datasets_scale[independPo]
datasets_nonscale_independ<-datasets_nonscale[independPo]
input_files<-names(datasets_nonscale)[independPo]

model7score.independ<-list()
for(i in 1:length(model.coef1)){
  modelGlobalScore.i<-c()
  for(j in 1:length(input_files)){
    modelGlobalScore.i<-c(modelGlobalScore.i,
                          modelScore(model_files[i],input_files[j]))
  }
  model7score.independ[[names(model.coef1)[i]]]<-modelGlobalScore.i
  
}
model7score.independ[["ssgsea"]]<-score.test.independ

#ROC3_5(scale(score),clinicalMatrixindepend,1,2)
ppp<-c()
for(i in 1:length(model7score.independ)){
  score.i<-model7score.independ[[i]]
  phe.i<-clinicalMatrixIndepend[match(names(score.i),rownames(clinicalMatrixIndepend)),]
  ppp<-rbind(ppp,ROC3_5_CI(scale(score.i),phe.i,1,2))
}
rownames(ppp)<-names(model7score.independ)
colnames(ppp)<-c("3year","5year")
allROC[["ssgsea_cox"]][["independ"]]<-ppp
#####test总######
testallPo<-match(c(names(datasets_nonscale)[testPo],
                   names(datasets_nonscale)[independPo]),
                 names(datasets_nonscale))
datasets_scale_testAll<-datasets_scale[testallPo]
datasets_nonscale_testAll<-datasets_nonscale[testallPo]
input_files<-names(datasets_nonscale)[testallPo]

model7score.testAll<-list()
for(i in 1:length(model.coef1.cox)){
  modelGlobalScore.i<-c()
  for(j in 1:length(input_files)){
    modelGlobalScore.i<-c(modelGlobalScore.i,
                          modelScore(model_files[i],input_files[j]))
  }
  model7score.testAll[[names(model.coef1.cox)[i]]]<-modelGlobalScore.i
  
}
model7score.testAll[["ssgsea"]]<-c(score.test,score.test.independ)

#ROC3_5(scale(score),clinicalMatrixtestAll,1,2)
ppp<-c()
for(i in 1:length(model7score.testAll)){
  score.i<-model7score.testAll[[i]]
  phe.i<-clinicalMatrixTestAll[match(names(score.i),rownames(clinicalMatrixTestAll)),]
  ppp<-rbind(ppp,ROC3_5(scale(score.i),phe.i,1,2))
}
rownames(ppp)<-names(model7score.testAll)
colnames(ppp)<-c("3year","5year")
allROC[["ssgsea_cox"]][["testAll"]]<-ppp

#####all samples######
allPo<-match(c(names(datasets_nonscale)[trainPo],
               names(datasets_nonscale)[testPo],
               names(datasets_nonscale)[independPo]),
             names(datasets_nonscale))
datasets_scale_All<-datasets_scale[allPo]
datasets_nonscale_All<-datasets_nonscale[allPo]
input_files<-names(datasets_nonscale)[allPo]

model15score.All<-list()
for(i in 1:length(model.coef1)){
  modelGlobalScore.i<-c()
  for(j in 1:length(input_files)){
    modelGlobalScore.i<-c(modelGlobalScore.i,
                          modelScore(model_files[i],input_files[j]))
  }
  model7score.All[[names(model.coef1.cox)[i]]]<-modelGlobalScore.i
  
}
model7score.All[["ssgsea"]]<-c(score,score.test,score.test.independ)

#ROC3_5(scale(score),clinicalMatrixAll,1,2)
ppp<-c()
for(i in 1:length(model7score.All)){
  score.i<-model7score.All[[i]]
  phe.i<-clinicalMatrixAll[match(names(score.i),rownames(clinicalMatrixAll)),]
  ppp<-rbind(ppp,ROC3_5(scale(score.i),phe.i,1,2))
}
rownames(ppp)<-names(model7score.All)
colnames(ppp)<-c("3year","5year")
allROC[["ssgsea_cox"]][["All"]]<-ppp
#####all samples-ci#####
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

######3,5year-meta#####
y3_matrix_meta<-c()
for(i in 1:nrow(ROC[[1]])){
  matrix_premeta<-c()
  for(j in 1:length(ROC)){
    matrix_premeta<-rbind(matrix_premeta,ROC[[j]][i,1:2])
  }
  #求C-statics的meta结果
  AUC<-meta.summaries(matrix_premeta[,1],matrix_premeta[,2],method="fixed",logscale=F)
  summary_AUC<-c()
  if(AUC$het[3]>0.05){
    summary_AUC<-c(AUC$summary,AUC$se.summary)
  }else{
    AUC<-meta.summaries(matrix_premeta[,1],matrix_premeta[,2],method="random",logscale=F)
    summary_AUC<-c(AUC$summary,AUC$se.summary)
  }
  y3_matrix_meta<-rbind(y3_matrix_meta,summary_AUC)
}

y5_matrix_meta<-c()
for(i in 1:nrow(ROC[[1]])){
  matrix_premeta<-c()
  for(j in c(1:length(ROC))[-c(3,11)]){
    matrix_premeta<-rbind(matrix_premeta,ROC[[j]][i,3:4])
  }
  #求C-statics的meta结果
  AUC<-meta.summaries(matrix_premeta[,1],matrix_premeta[,2],method="fixed",logscale=F)
  summary_AUC<-c()
  if(AUC$het[3]>0.05){
    summary_AUC<-c(AUC$summary,AUC$se.summary)
  }else{
    AUC<-meta.summaries(matrix_premeta[,1],matrix_premeta[,2],method="random",logscale=F)
    summary_AUC<-c(AUC$summary,AUC$se.summary)
  }
  y5_matrix_meta<-rbind(y5_matrix_meta,summary_AUC)
}


Multi_matrix_meta<-cbind(y3_matrix_meta,y5_matrix_meta)

rownames(Multi_matrix_meta)<-rownames(ROC[[1]])
colnames(Multi_matrix_meta)<-c("aucOf3","se.aucof3","aucOf5","se.aucof5")

#graph
Data<-Multi_matrix_meta
Data<-as.data.frame(Data[order(Data[,1],decreasing = T),])
Data$model<-factor(rownames(Data),levels =rownames(Data) )
#Data$model<-rownames(Data)
# ggplot(Data, aes(x=model, y=aucOf3)) + 
#   geom_bar(stat ="identity",width = 0.65,fill="#BE0005") +
#   geom_errorbar(aes(ymin=aucOf3-1.96*se.aucof3, ymax=aucOf3+1.96*se.aucof3), width=.1) +
#   #geom_errorbar(aes(ymin=aucOf3-se.aucof3, ymax=aucOf3+se.aucof3), width=.1)
#   labs(y="Survival ROC")+
#   scale_y_continuous(expand = c(0,0),breaks = seq(0,0.65,0.05))+
#   coord_cartesian(ylim = c(0.4, 0.7))+
#   theme_bw()+
#   theme(panel.border = element_blank(),panel.grid=element_blank(),
#         axis.line = element_line(colour = "black"),
#         axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1,color = "black",
#                                    size=14),
#         axis.text.y=element_text(color="black",size=14),
#         plot.title = element_text(colour = "black", face = "bold", 
#                                   size = 14, vjust = 1,hjust = 0.5),
#         axis.title.x = element_blank(),
#         axis.title.y=element_text(color="black",size=16),
#         legend.position ="none")
# ggsave("allmodels_average_ROC.pdf",width = 7,height=5)


meta.i<-metagen(aucOf3,se.aucof3,sm="AUC",data = Data,
                studlab = rownames(Data))
pdf("allmodels'forest_metaROC.pdf",width=7,height = 7)
forest(meta.i, leftcols = c("studlab"),leftlabs = c("Study Element"),
       rightcols = c("effect.ci"),fontsize = 13,just.studlab = "center",
       overall = F,col.square = "#363886",lwd=2,xlim = c(0.4,0.7),
       at=c(0.4,0.5,0.6,0.7),plotwidth="3.5cm",squaresize = 0.4,
       weight.study = "same",col.inside="black",ref = 0.5)
dev.off()

heatmapData<-c()
for(i in 1:length(ROC)){
  heatmapData<-cbind(heatmapData,ROC[[i]][,1])
}
colnames(heatmapData)<-names(ROC)

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
  geom_bar(stat ="identity",width = 0.8,fill="#BE0005") +#"#EFC000"
  geom_errorbar(aes(ymin=aucOf3-1.96*se.aucof3, ymax=aucOf3+1.96*se.aucof3), width=.1) +
  #geom_errorbar(aes(ymin=aucOf3-se.aucof3, ymax=aucOf3+se.aucof3), width=.1)
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


# 
# #AUCname<-rownames(Data)
# meta.i<-metagen(aucOf3,se.aucof3,sm="AUC",data = Data,
#                 studlab = rownames(Data))
# forest(meta.i, leftcols = c("studlab"),leftlabs = c("Study Element"),
#        rightcols = c("effect.ci"),fontsize = 13,just.studlab = "center",
#        overall = F,col.square = "#363886",lwd=2,xlim = c(0.4,0.7),
#        at=c(0.4,0.5,0.6,0.7),plotwidth="3.5cm",squaresize = 0.4,
#        weight.study = "same",col.inside="black",ref = 0.5)


######3,5year-average#####
y3_matrix_average<-c()
for(i in 1:nrow(ROC[[1]])){
  matrix_premeta<-c()
  for(j in 1:length(ROC)){
    matrix_premeta<-rbind(matrix_premeta,ROC[[j]][i,1:2])
  }
  #AUC<-apply(matrix_premeta,2,mean)
  AUC<-c(mean(matrix_premeta[,1]),sd(matrix_premeta[,1])/sqrt(15))
  y3_matrix_average<-rbind(y3_matrix_average,AUC)
}

y5_matrix_average<-c()
for(i in 1:nrow(ROC[[1]])){
  matrix_premeta<-c()
  for(j in c(1:length(ROC))[-c(3,11)]){
    matrix_premeta<-rbind(matrix_premeta,ROC[[j]][i,3:4])
  }
  #AUC<-apply(matrix_premeta,2,mean)
  AUC<-c(mean(matrix_premeta[,1]),sd(matrix_premeta[,1])/sqrt(15))
  y5_matrix_average<-rbind(y5_matrix_average,AUC)
}


Multi_matrix_average<-cbind(y3_matrix_average,y5_matrix_average)

rownames(Multi_matrix_average)<-rownames(ROC[[1]])
colnames(Multi_matrix_average)<-c("aucOf3","se.aucof3","aucOf5","se.aucof5")



#graph
Data<-Multi_matrix_average
Data<-as.data.frame(Data[order(Data[,1],decreasing = T),])
Data$model<-factor(rownames(Data),levels =rownames(Data) )
#Data$model<-rownames(Data)
ggplot(Data, aes(x=model, y=aucOf3)) + 
  geom_bar(stat ="identity",width = 0.65,fill="#BE0005") +
  geom_errorbar(aes(ymin=aucOf3-1.96*se.aucof3, ymax=aucOf3+1.96*se.aucof3), width=.1) +
  #geom_errorbar(aes(ymin=aucOf3-se.aucof3, ymax=aucOf3+se.aucof3), width=.1)
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
ggsave("allmodels_average_ROC.pdf",width = 7,height=5)

# 
# meta.i<-metagen(aucOf3,se.aucof3,sm="AUC",data = Data,
#                 studlab = rownames(Data))
# forest(meta.i, leftcols = c("studlab"),leftlabs = c("Study Element"),
#        rightcols = c("effect.ci"),fontsize = 13,just.studlab = "center",
#        overall = F,col.square = "#363886",lwd=2,xlim = c(0.4,0.7),
#        at=c(0.4,0.5,0.6,0.7),plotwidth="3.5cm",squaresize = 0.4,
#        weight.study = "same",col.inside="black",ref = 0.5)

heatmapData<-c()
for(i in 1:nrow(ROC[[1]])){
  heatmapData<-cbind(heatmapData,ROC[[i]][,1])
}
colnames(heatmapData)<-zz

heatmapData<-heatmapData[match(rownames(Data),rownames(heatmapData)),]
bk<-c(seq(0.3,0.49,length.out=19),0.5,seq(0.51,0.81,length.out=30))
color<-c(colorRampPalette(c("#7A7ABC","white"))(19),"white",
         colorRampPalette(c("white","firebrick3"))(30),"firebrick3")
pdf("allmodels_average_ROC-heatmap.pdf",width=6,height=5)
pheatmap(heatmapData,cluster_cols = F,cluster_rows = F,
         color=color,breaks=bk,border_color="white",
         cellwidth = 15,cellheight = 15)
dev.off()



Data<-Multi_matrix_average[c(1,2,3,5,6,7,15),]
Data<-as.data.frame(Data[order(Data[,1],decreasing = T),])
Data$model<-factor(rownames(Data),levels =rownames(Data) )
ggplot(Data, aes(x=model, y=aucOf3)) + 
  geom_bar(stat ="identity",width = 0.8,fill="#BE0005") +#"#EFC000"
  geom_errorbar(aes(ymin=aucOf3-1.96*se.aucof3, ymax=aucOf3+1.96*se.aucof3), width=.1) +
  #geom_errorbar(aes(ymin=aucOf3-se.aucof3, ymax=aucOf3+se.aucof3), width=.1)
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
ggsave("coxmodels_average_ROC.pdf",width = 7,height=7)


# 
# #AUCname<-rownames(Data)
# meta.i<-metagen(aucOf3,se.aucof3,sm="AUC",data = Data,
#                 studlab = rownames(Data))
# forest(meta.i, leftcols = c("studlab"),leftlabs = c("Study Element"),
#        rightcols = c("effect.ci"),fontsize = 13,just.studlab = "center",
#        overall = F,col.square = "#363886",lwd=2,xlim = c(0.4,0.7),
#        at=c(0.4,0.5,0.6,0.7),plotwidth="3.5cm",squaresize = 0.4,
#        weight.study = "same",col.inside="black",ref = 0.5)
