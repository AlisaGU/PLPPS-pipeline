stage_grade_age<-function(phe,gradeStatus,stageColname,gradeColname,ageColname){
  #phe是挑选样本之前的临床信息表
  #gradeStatus表示grade这一列是否有非NA值
  if(gradeStatus==0){
    #表示全为NA
    phe<-phe[which(phe$histological_type=="ser" & phe$sample_type=="tumor" & phe$summarystage=="late"),]
  }else if(gradeStatus==1){
    phe<-phe[which(phe$histological_type=="ser" & phe$sample_type=="tumor" & phe$summarystage=="late" & phe$grade>=2),]
  }else{
    phe<-phe
  }
  #之后的phe是晚期，HG样本
  stage.col<-phe[,which(colnames(phe)==stageColname)]
  grade.col<-phe[,which(colnames(phe)==gradeColname)]
  age.col<-phe[,which(colnames(phe)==ageColname)]
  if(any(!is.na(stage.col))){
    num<-length(which(stage.col>3))
    ratio<-format(num/length(stage.col)*100,digits=3)
    stage<-paste(num,"(",ratio,")",sep= "")
  }else{stage<-NA}
  if(any(!is.na(grade.col))){
    num<-length(which(grade.col>=3))
    ratio<-format(num/length(grade.col)*100,digits=3)
    grade<-paste(num,"(",ratio,")",sep= "")
  }else{grade<-NA}
  if(any(!is.na(age.col))){
    age<-age.col[!is.na(age.col)]
    meanAge<-mean(age)
    SD<-sd(age)
    confidenceInterval<-1.96*SD/sqrt(length(age))
    low<-format(meanAge-confidenceInterval,digits=3)
    high<-format(meanAge+confidenceInterval,digits=3)
    age<-paste(format(meanAge,digits=3),"(",low,"-",high,")",sep= "")
  }else{age<-NA}
  return(c(nrow(phe),stage,grade,age))#第一个为晚期、HG样本的数目
}