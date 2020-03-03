# Filename: PErrC.r
# Title: Determine prediction error curve for a test set.
#
# Author: Christoph Bernau
# Email: <bernau@ibe.med.uni-muenchen.de>
# Date of creation: April 27, 2012
#
# Brief description:
#   Returns vector of prediction errors corresponding to the timegrid in 
#	prob.breslow.
#
# Arguments:
#   -prob.breslow: breslow-object returned by function "getSurvProbs"
#   -newy: survival response for new observations
#	 -y: training survival response (used for calculating censoring 
#	probabilities)
#
# Value:
#   vector of prediction errors 
#
###############################################################################



setMethod("evaluate", signature(object = "PErrC"), function(object, 
                                                            prob.breslow, newy, add = list(),...) {
  if(is.null(add$ytrain)) stop('Survival Response for training data is 
                               missing (argument "add" in evaluate-PErrC).')	
  ytrain<-add$ytrain
  yc<-rbind(ytrain,newy)
  timegrid <- prob.breslow@time
  # get breslow estimator for censoring on full data set 
  #(see peperr pmpec)
  censdata <- data.frame(time = yc[, 1], status = 
                           ifelse(yc[, 2] == 0, 1, 0))
  cens.fit <- survival::survfit(Surv(time, status) 
                                ~ 1, data = censdata)
  # evaluate breslow estimator on timegrid
  eval.cens.prob <- summary(cens.fit, times = timegrid)$surv
  
  if (length(eval.cens.prob) != length(timegrid)) 
    stop("Censoring Probability cannot be estimated for some 
         timepoints  in timegrid. Probably, the last ", 
         length(timegrid) - length(eval.cens.prob),
         " are too late.")
  
  ### needed for easier computation lateron
  sort.time <- sort(yc[, 1])
  sortsurv <- c(1, summary(cens.fit, times = sort.time)$surv)
  sort.time <- c(-999, sort.time)
  sortsurv <- c(1, sortsurv)
  
  invcensprob <- unlist(lapply(newy[, 1], function(actual.time) 
    sortsurv[which(sort.time >= 
                     actual.time)[1] - 1]))
  
  time <- newy[, 1]
  status <- newy[, 2]
  
  ## statusmatrix for new obs on timegrid
  status.mat <- matrix(timegrid, length(time), length(timegrid), 
                       byrow = TRUE) < time
  ## weightatrix for new obs on timegrid
  weight.status.mat <- matrix(timegrid, length(time), 
                              length(timegrid), byrow = TRUE) < time
  
  ### weights as in err(bo,cv) in binder 2010
  weightmat <- t(t(weight.status.mat)/eval.cens.prob) +
    (1 - weight.status.mat) * matrix(status != 0, 
                                     length(status), length(timegrid)) * 
    matrix(1/invcensprob, 
           length(status), length(timegrid))
  
  ### compute err(bo,cv)
  return( apply(weightmat * (status.mat - prob.breslow@curves)^2, 
                2, mean) )
})

### aggregate method
setMethod("aggregateMeasure", signature(object = "PErrC"), 
          function(object, foldmeasures) {
            return(mean(unlist(foldmeasures), na.rm = TRUE))
          })
