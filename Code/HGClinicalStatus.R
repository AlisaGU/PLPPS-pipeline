HGClinicalStatus<-function(data,gradeStatus){
  idx.tumorstage<-structure(rep(NA,length(data)),names=names(data))
  idx.tumorgrade<-structure(rep(NA,length(data)),names=names(data))
  idx.age<-structure(rep(NA,length(data)),names=names(data))
  for(i in 1:length(data)){
    X<-data[[i]]
    if(gradeStatus[i]==0){
      idx.tumorstage[i] <- (sum(!is.na(X$tumorstage[X$summarystage=="late"])) > 0)& (length(unique(X$tumorstage[X$summarystage=="late"])) > 1)
      idx.tumorgrade[i] <- (sum(!is.na(X$grade[X$summarystage=="late"])) > 0) & (length(unique(X$grade[X$summarystage=="late"])) > 1)
      idx.age [i]<- (sum(!is.na(X$age_at_initial_pathologic_diagnosis[X$summarystage=="late"])) > 0)& (length(unique(X$age_at_initial_pathologic_diagnosis[X$summarystage=="late"])) > 1)
    }else if(gradeStatus[i]==1){
      idx.tumorstage[i] <- (sum(!is.na(X$tumorstage[X$summarystage=="late"&X$grade>=2])) > 0) & (length(unique(X$tumorstage[X$summarystage=="late"&X$grade>=2])) > 1)
      idx.tumorgrade[i] <- (sum(!is.na(X$grade[X$summarystage=="late"&X$grade>=2])) > 0) & (length(unique(X$grade[X$summarystage=="late"&X$grade>=2])) > 1)
      idx.age [i]<-(sum(!is.na(X$age_at_initial_pathologic_diagnosis[X$summarystage=="late"&X$grade>=2])) > 0)& (length(unique(X$age_at_initial_pathologic_diagnosis[X$summarystage=="late"&X$grade>=2])) > 1)
    }else if(gradeStatus[i]==2){
      idx.tumorstage[i] <- T
      idx.tumorgrade[i] <- T
      idx.age [i]<-T
    }
  }
  return(cbind(idx.tumorstage,idx.tumorgrade,idx.age))
}
