# Filename: cvPLogL.r
# Title: Determine partial log likelihood for a testset.
#
# Author: Christoph Bernau
# Email: <bernau@ibe.med.uni-muenchen.de>
# Date of creation: April 27, 2012
#
# Brief description:
#   Returns partial loglikelihood for the testset.
#
# Arguments
#   -pred.all: linear predictor for all observations
#   -yc: complete survival response (class Surv) 
#   -learnind:vector indicating which observations are to be used for training
#	the survival model
#
# Value:
#   cross-validated partial log likelihood for the fold specified via learnind
#	(class numeric)
#
###############################################################################



setMethod("evaluate", signature(object = "CvPLogL"), function(object, pred.test, 
                                                              newy,add = list(),...) {
  
  if(is.null(add$ytrain)) stop('Survival Response for training data is 
                               missing (argument "add" in ev-CvPLogL).')	
  if(is.null(add$lptrain)) stop('Linear Predictor for training data is 
                                missing (argument "add" in ev-CvPLogL).')
  
  pred.all<-c(add$lptrain@lp,pred.test)
  yc<-rbind(add$ytrain,newy)
  
  
  ## Verweij Cross-validated log partial likelihood
  full.loglik <- .logPLik(pred.all, yc)
  fold.loglik <- .logPLik(add$lptrain@lp, add$ytrain)
  cvi <- full.loglik - fold.loglik
  cvi
})

setMethod("aggregateMeasure", signature(object = "CvPLogL"), function(object, 
                                                                      foldmeasures) {
  mean(unlist(foldmeasures))
})

.logPLik <- function(linear.predictor, y) {
  ## adapted from GPL code from the CoxBoost CRAN package
  efron.weightmat <- function(time, status) {
    n <- length(time)
    uncens <- which(status == 1)
    weightmat <- matrix(0, n, length(uncens))
    rept <- rep(0, n)
    for (i in 1:n) rept[i] <- sum(time[i:n] == time[i] & status[i:n] == 1)
    for (i in 1:length(uncens)) {
      weightmat[time >= time[uncens][i] |
                  (status != 0 & status != 1), i] <- 1
      tie <- time == time[uncens[i]] & status == 1
      di <- max(rept[tie])
      weightmat[tie, i] <- weightmat[tie, i] - (di - rept[uncens[i]])/di
    }
    return(weightmat)
  }
  newtime <- as.numeric(y)[1:nrow(y)]
  newstatus <- as.numeric(y)[-1:-nrow(y)]
  uncens <- which(newstatus == 1)
  weightmat <- efron.weightmat(newtime, newstatus)
  logplik <- sum(linear.predictor[uncens] - 
                   log(apply(weightmat * exp(linear.predictor), 2, sum)))
  return(logplik)
}

