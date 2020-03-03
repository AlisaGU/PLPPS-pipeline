# Filename: gbmbasehaz.r
# Author: C. Bernau
# Email: <bernau@ibe.med.uni-muenchen.de>
# Date of creation: 16.5.2012
#
# Brief description:
#   predict baselinehazard on a timegrid for ModelLearnedobjects
#
# Arguments
#   -survout: output-object of a *surv-function (class ModelLearned)
#   -timegrid: grid of timepoints for which the probabilities are to 
#	be predicted
#
# Value:
#   a vector of predicted baseline-hazards for the specified timegrid
#
###############################################################################

gbmbasehaz<-function(survout,timegrid){
  require(gbm)
  learnind<-survout@learnind
  oldy <- survout@y
  oldlp <- survout@linear.predictor@lp
  if(length(oldlp)==0) stop("Model does not provide a linear predictor.")
  base<-basehaz.gbm(t=oldy[,1],delta=oldy[,2],t.eval=timegrid,f.x=oldlp,
                    cumulative=TRUE,smooth=FALSE)
  ##t=0 should have base=0, not NA.  Note all.equal has tolerance for 
  ###numeric imprecision:
  if(all.equal(timegrid[1],0.0)==TRUE & is.na(base[1]))
    base[1] <- 0
  ##the 'base' survival function
  survbase<-exp(-base)
  survbase
}
