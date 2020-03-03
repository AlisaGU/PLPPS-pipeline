Train_meta_summary1<-Train_meta_summary
colnames(Train_meta_summary1)<-c("Q","df.Q","p.heterogeneity","I2",
                                 "HR","upper","lower","P",
                                 "HR","upper","lower","P")
gene_premodel<-c()

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


geneNamePreModel<-rownames(gene_enrich)


h<-enrichGO(geneNamePreModel,
            OrgDb="org.Hs.eg.db",
            keytype = "SYMBOL",
            ont="BP",
            pvalueCutoff = 1,
            pAdjustMethod = "fdr",
            minGSSize = 500,
            maxGSSize = 2000)
enrich_meta<-h[h[,7]<0.001,]

enrich_gene<-unique(unlist(strsplit(enrich_meta$geneID,"/")))


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


