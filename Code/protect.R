protect<-function(data){#data存储meta统计量
  protectGeneInfor<-c()#存储连续情况下置信区间小于1的基因名称
  for(i in 1:nrow(data)){
    x<-unlist(data[i,])
    if(x[3]<0.1&x[10]<1)#即拒绝零假设，认为具有异质性
    {
      protectGeneInfor<-rbind(protectGeneInfor,c(x[9:12],2))#2为随机效应模型
      rownames(protectGeneInfor)[nrow(protectGeneInfor)]<-rownames(data)[i]
    }
    if(x[3]>=0.1&x[6]<1)#接受零假设，认为具有同质性
    {
      protectGeneInfor<-rbind(protectGeneInfor,c(x[5:8],1))#1为固定效应模型
      rownames(protectGeneInfor)[nrow(protectGeneInfor)]<-rownames(data)[i]
    }
  }
  colnames(protectGeneInfor)<-c("HR","upper","lower","p value","model")
  protectGeneInfor<-protectGeneInfor[which(protectGeneInfor[,4]<0.01),]
  protectGene<-rownames(protectGeneInfor)
  return(list(protectGeneInfor=protectGeneInfor,
              protectGene=protectGene))
}