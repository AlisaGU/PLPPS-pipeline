getSurvivalMatrix1<-function(clinicalData,type,exceptSet,DataName=NULL){
  clinicalMatrix<-c()
  if(type=="global"){
    for(i in 1:length(clinicalData)){
      if(names(clinicalData)[i]!=exceptSet){
        preRow<-rownames(clinicalMatrix)
        clinicalMatrix<-rbind(clinicalMatrix,cbind(clinicalData[[i]]$days_to_death,
                                                   clinicalData[[i]]$vital_status))
        rownames(clinicalMatrix)<-c(preRow,rownames(clinicalData[[i]]))
      }else{
        preRow<-rownames(clinicalMatrix)
        clinicalMatrix<-rbind(clinicalMatrix,cbind(clinicalData[[i]]$donor_survival_time,
                                                   clinicalData[[i]]$donor_vital_status))
        rownames(clinicalMatrix)<-c(preRow,rownames(clinicalData[[i]]))
      }
    }
  }else if(type=="seperate"){
    if(DataName!=exceptSet){
      preRow<-rownames(clinicalMatrix)
      clinicalMatrix<-rbind(clinicalMatrix,cbind(clinicalData$days_to_death,
                                                 clinicalData$vital_status))
      rownames(clinicalMatrix)<-c(preRow,rownames(clinicalData))
    }else{
      preRow<-rownames(clinicalMatrix)
      clinicalMatrix<-rbind(clinicalMatrix,cbind(clinicalData$donor_survival_time,
                                                 clinicalData$donor_vital_status))
      rownames(clinicalMatrix)<-c(preRow,rownames(clinicalData))
    }
  }
  colnames(clinicalMatrix)<-c("time","status")
  for(i in 1:nrow(clinicalMatrix)){
    if(clinicalMatrix[i,2]=="1")
      clinicalMatrix[i,2]<-"Dead"
    if(clinicalMatrix[i,2]=="0")
      clinicalMatrix[i,2]<-"Alive"
  }
  return(clinicalMatrix)
}