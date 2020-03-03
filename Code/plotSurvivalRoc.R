plotSurvivalRoc<-function(SR,title){
  plot(SR$FP, SR$TP, type="l", xlim=c(0,1), ylim=c(0,1),   
       xlab=paste( "FP", "\n", "AUC = ",round(SR$AUC,3)), 
       ylab="TP",main=title)
  abline(0,1)
}