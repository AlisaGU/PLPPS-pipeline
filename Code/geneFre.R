geneFre<-function(Data){#Data存储训练集基因的lnHR,HR,se(lnHR)
  genes_yuanshi <- c()
  for (i in 1:length(Data)) {
    gene.meta.i<-Data[[i]]; 
    genes_yuanshi <- c(genes_yuanshi,rownames(gene.meta.i))
  }
  genes_yuanshi <- table(genes_yuanshi);
  genes_frequency <- genes_yuanshi[genes_yuanshi>2] 
  
  gene.metaL2<-list()#存储满足频率条件的基因的lnHR，HR,se(lnHR)
  for(i in 1:length(names(genes_frequency))){
    gene.i<-names(genes_frequency)[i]
    gene.metaL2.i<-c()
    for(j in 1:length(Data)){
      if(is.element(gene.i,rownames(Data[[j]]))){
        position<-which(rownames(Data[[j]])==gene.i)
        gene.metaL2.i<-rbind(gene.metaL2.i,Data[[j]][position,])
        rownames(gene.metaL2.i)[nrow(gene.metaL2.i)]<-names(Data)[j]
      }
    }
    gene.metaL2[[i]]<-gene.metaL2.i
  }
  names(gene.metaL2)<-names(genes_frequency)
  return(gene.metaL2)
}