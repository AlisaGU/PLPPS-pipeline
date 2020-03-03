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