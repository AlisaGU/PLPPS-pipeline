protect<-function(data){#data�洢metaͳ����
  protectGeneInfor<-c()#�洢�����������������С��1�Ļ�������
  for(i in 1:nrow(data)){
    x<-unlist(data[i,])
    if(x[3]<0.1&x[10]<1)#���ܾ�����裬��Ϊ����������
    {
      protectGeneInfor<-rbind(protectGeneInfor,c(x[9:12],2))#2Ϊ���ЧӦģ��
      rownames(protectGeneInfor)[nrow(protectGeneInfor)]<-rownames(data)[i]
    }
    if(x[3]>=0.1&x[6]<1)#��������裬��Ϊ����ͬ����
    {
      protectGeneInfor<-rbind(protectGeneInfor,c(x[5:8],1))#1Ϊ�̶�ЧӦģ��
      rownames(protectGeneInfor)[nrow(protectGeneInfor)]<-rownames(data)[i]
    }
  }
  colnames(protectGeneInfor)<-c("HR","upper","lower","p value","model")
  protectGeneInfor<-protectGeneInfor[which(protectGeneInfor[,4]<0.01),]
  protectGene<-rownames(protectGeneInfor)
  return(list(protectGeneInfor=protectGeneInfor,
              protectGene=protectGene))
}