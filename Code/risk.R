risk<-function(data){#data存储meta统计量
  riskGeneInfor<-c()#存储连续情况下置信区间大于1的基因名称
  for(i in 1:nrow(data)){
    x<-unlist(data[i,])
    if(x[3]<0.1&x[11]>1)#即拒绝零假设，认为具有异质性
    {
      riskGeneInfor<-rbind(riskGeneInfor,c(x[9:12],2))#2为随机效应模型
      rownames(riskGeneInfor)[nrow(riskGeneInfor)]<-rownames(data)[i]
    }
    if(x[3]>=0.1&x[7]>1)#接受零假设，认为具有同质性
    {
      riskGeneInfor<-rbind(riskGeneInfor,c(x[5:8],1))#1为固定效应模型
      rownames(riskGeneInfor)[nrow(riskGeneInfor)]<-rownames(data)[i]
    }
  }
  colnames(riskGeneInfor)<-c("HR","upper","lower","p value","model")
  riskGeneInfor<-riskGeneInfor[which(riskGeneInfor[,4]<0.01),]
  riskGene<-rownames(riskGeneInfor)
  return(list(riskGeneInfor=riskGeneInfor,
              riskGene=riskGene))
}