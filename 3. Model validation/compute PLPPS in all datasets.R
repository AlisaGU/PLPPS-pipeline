source('/code/double-weighted_ssGSEA.R')
source('/code/ssgseaScore.R')
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
                            