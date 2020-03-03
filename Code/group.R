group<-function(x,type,median){
  #numChara是一个向量，分别为（min,1st Qu,median, 3rd Qu,max）
  Group<-""
  if(type=="summary"){
    if(x>=numChara[1]&x<=numChara[2])
      Group<-"one"
    if(x>=numChara[2]&x<=numChara[3])
      Group<-"two"
    if(x>=numChara[3]&x<=numChara[4])
      Group<-"three"
    if(x>=numChara[4]&x<=numChara[5])
      Group<-"four"
  }else{
    if(x>=median)
      Group<-"negative"#High Risk
    if(x<median)
      Group<-"positive"
  }
  return(Group)
}