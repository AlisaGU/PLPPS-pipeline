ssgseaScore<-function(Data,w1,gProtect,gRisk,wProtect,wRisk){
  #将ssgsea模型应用到Data列表的每一个元素。
  ssgsea_out<-list()
  for(i in 1:length(Data))
  {
    protect<-unlist(Project.to.GeneSet(data.array=as.matrix(Data[[i]]),
                                    weight=w1,gene.set=gProtect,w.gset=wProtect))
    risk<-unlist(Project.to.GeneSet(data.array=as.matrix(Data[[i]]),
                                   weight=w1,gene.set=gRisk,w.gset=wRisk))
    sScore<-risk-protect
    scaledsScore<-sScore/nrow(Data[[i]])#除以数据集的基因总数，进行标准化
    ssgsea_out[[i]]<-rbind(protect,risk,sScore,scaledsScore)
    colnames(ssgsea_out[[i]])<-colnames(Data[[i]])
    rownames(ssgsea_out[[i]])<-c("Protect Score","Risk Score",
                                 "nonscaled Global Score","scaled Global Score")
  }
  names(ssgsea_out)<-names(Data)
  return(ssgsea_out)
}