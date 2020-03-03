setwd("")#data path
source('/code/HGClinicalStatus.R')
source('/code/coxAnalysis.R')
library(GEOquery)
library(survival) 
load("clinical.Rdata")
load("HGSC datasets with OS.rda")
gradeStatus<-c(0,rep(1,5),0,rep(1,8),2)
gradeStatus[14]<-2
result<-HGClinicalStatus(clinical,gradeStatus)
idx.tumorstage<-result[,1]
idx.tumorgrade<-result[,2]
idx.age<-result[,3]
expNa<-c(rep(0,13),1,1)
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

load("gene_cox.Rdata")
threshold<-0.01
actual_genenum<-matrix(sapply(gene.cox,function(x){length(which(x[,1]<threshold))}),nrow=1)
rownames(actual_genenum)<-"genes related to survival"
colnames(actual_genenum)<-names(gene.cox)
diff.coxgene<-lapply(gene.cox,function(x){rownames(x)[which(x[,1]<threshold)]})

source('/code/lower2symmetric.R')
source('/codes/signatureSimility.R')
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

source('/code/HGClinicalStatus.R')
source('/code/coxAnalysis.R')
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

