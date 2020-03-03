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

protectModelGene<-as.character(unlist(read.table("protectModelGene.txt",as.is=T)))
riskModelGene<-as.character(unlist(read.table("riskModelGene.txt",as.is=T)))
protectGeneInfor<-read.csv("protectGeneInfor.csv",as.is=T,header = T,row.names = 1)
riskGeneInfor<-read.csv("riskGeneInfor.csv",as.is=T,header = T,row.names = 1)
protectWeight<-abs(log(protectGeneInfor[protectModelGene,1]))
riskWeight<-abs(log(riskGeneInfor[riskModelGene,1]))
