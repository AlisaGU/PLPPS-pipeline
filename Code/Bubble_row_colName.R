Bubble_row_colName<-function(data,pos){#返回上、下三角矩阵每一个值对应的横纵坐标
  Name<-c()
  if(pos=="lower"){
    lowerMa<-lower.tri(data)
    for(i in 2:nrow(data)){
      for(j in 1:(i-1) )
        Name<-rbind(Name,c(rownames(data)[i],colnames(data)[j]))
    }
  }
  if(pos=="upper"){
    upperMa<-upper.tri(data)
    for(i in 1:(nrow(data)-1)){
      for(j in (i+1):ncol(data)){
        Name<-rbind(Name,c(rownames(data)[i],colnames(data)[j]))
      }
    }
  }
  if(pos=="all"){#by col
  for(i in 1:nrow(data))
    for(j in 1:ncol(data))
      Name<-rbind(Name,c(rownames(data)[i],colnames(data)[j]))
  }
  o<-data.frame(x=vector(length=nrow(data)^2),
                y=vector(length=nrow(data)^2))
  o$x<-factor(Name[,1],levels=rownames(data))
  o$y<-factor(Name[,2],levels=rownames(data))
  return(o)
}