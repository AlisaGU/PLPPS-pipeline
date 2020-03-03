
source('/code/geneFre.R')
load("datasets_scale.Rdata")
load("gene_meta.Rdata")
trainPo<-which(sapply(datasets_scale,ncol)>100)[-8]
trainSet<-gene.meta[trainPo]
gene.metaTrain<-geneFre(trainSet)




source('/code/meta_summary.R')
library("meta")
library(clusterProfiler)
library(org.Hs.eg.db)
load("gene_metaTrain.Rdata")
Train_meta_summary<-meta_summary(gene.metaTrain)




source('/code/risk.R')
source('/code/protect.R')
Train_meta_summary<-read.csv("Train_meta_summary.csv",
                             header = T,row.names = 1,as.is=T)
colnames(Train_meta_summary)<-NULL
a<-risk(Train_meta_summary)
riskGeneInfor<-a[[1]]
riskGene<-a[[2]]
b<-protect(Train_meta_summary)
protectGeneInfor<-b[[1]]
protectGene<-b[[2]]
