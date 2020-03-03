source('/code/Bubble.R')
source('/code/Bubble_row_colName.R')
source('/code/bubblePlot.R')
library(ggplot2)


colnames(tri_semsim_term)<-rownames(tri_semsim_term)
colnames(tri_semsim_term_p)<-rownames(tri_semsim_term_p)
bubblePlot(tri_semsim_term,tri_semsim_term_p,
           "cox sim based on Term permutation",
           "cox sim based on Term permutation.pdf")
            
            
load("diff_coxgene.Rdata")
semSimilirity<-read.csv("semSimilirity.csv",header=T,row.names = 1,as.is=T)
semSimilirity_p<-read.csv("semSimilirity_p.csv",header=T,row.names = 1,as.is=T)
colnames(semSimilirity)<-colnames(semSimilirity_p)<-rownames(semSimilirity)<-rownames(semSimilirity_p)<-names(diff.coxgene)
bubblePlot(semSimilirity,semSimilirity_p,
           "cox sim based on gene permutation",
           "cox sim based on gene permutation.pdf")
           
load("diff_coxgene.Rdata")
coxGeneSim<-read.csv("coxGeneSim.csv",header=T,row.names = 1,as.is=T)
coxGeneSimP<-read.csv("coxGeneSimP.csv",header=T,row.names = 1,as.is=T)
colnames(coxGeneSim)<-colnames(coxGeneSimP)<-rownames(coxGeneSim)<-rownames(coxGeneSimP)<-names(diff.coxgene)
bubblePlot(coxGeneSim,coxGeneSimP,
           "cox sim based on exp correlation",
           "cox sim based on exp correlation.pdf")
           
load("diff_coxgene.Rdata")
semSimilirity<-read.csv("semSimilirity.csv",header=T,row.names = 1,as.is=T)
semSimilirity_p<-read.csv("semSimilirity_p.csv",header=T,row.names = 1,as.is=T)
colnames(semSimilirity)<-colnames(semSimilirity_p)<-rownames(semSimilirity)<-rownames(semSimilirity_p)<-names(diff.coxgene)
coxGeneSim<-read.csv("coxGeneSim.csv",header=T,row.names = 1,as.is=T)
coxGeneSimP<-read.csv("coxGeneSimP.csv",header=T,row.names = 1,as.is=T)
colnames(coxGeneSim)<-colnames(coxGeneSimP)<-rownames(coxGeneSim)<-rownames(coxGeneSimP)<-names(diff.coxgene)


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
           
jaccardMatrix<-read.csv("jaccardMatrix.csv",header=T,row.names = 1,as.is=T)
jaccardMatrix_p<-read.csv("jaccardMatrix_p.csv",header=T,row.names = 1,as.is=T)
colnames(jaccardMatrix)<-rownames(jaccardMatrix)
colnames(jaccardMatrix_p)<-rownames(jaccardMatrix_p)
bubblePlot(jaccardMatrix,jaccardMatrix_p,
           "cox sim based on jaccard",
           "cox sim based on jaccard.pdf")
           
jaccardMatrix<-read.csv("jaccardMatrix.csv",header=T,row.names = 1,as.is=T)
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
        