# Filename: getBestParameters.r
# Title: Determine best tuning parameter.
#
# Author: Christoph Bernau
# Email: <bernau@ibe.med.uni-muenchen.de>
# Date of creation: April 26, 2012
#
# Brief description:
#   Retruns list consisting of the index of the best parameter and the 
#	performances of all parameters.
#
# Arguments:
#   -tuneres: return object of method "tune"
#   -res.ind: index of resampling iteration
#   -measure: evaluation metric to be used
#
# Value:
#   list of index of best parameter and performances of all parameters
#   
###############################################################################


getBestParameters <- function(tuneres, res.ind, measure = NULL, ...) {
  
  if (is.null(measure)) {
    ## default tuning is cross-validated partial log-likelihood:
    measurec <- new("CvPLogL")
  } else {
    measurec <- new(measure)
  }
  
  tuneresi <- tuneres@tuneres[[res.ind]]
  
  if (class(tuneresi[[1]]) == "EvaluateOutput") 
    measurec <- new(tuneresi[[1]]@measure)
  
  ll <- list(...)
  ll$measure <- measure
  ll$X <- tuneres@X[tuneres@LearningSets@learnmatrix[res.ind, ], ]
  
  perf <- numeric(length(tuneresi))
  # Xcurr<-tuneres@X[tuneres@LearningSets@learnmatrix[res.ind,],]
  for (k in 1:length(tuneresi)) {
    
    if (class(tuneresi[[k]]) == "LearnOut") {
      ll$object <- tuneresi[[k]]
      foldmeasures <- do.call(evaluate, args = ll)@result  
      #tuneresi[[k]],measure=measure,X=Xcurr)@result
    } else foldmeasures <- tuneresi[[k]]@result
    perf[k] <- aggregateMeasure(measurec, foldmeasures)
  }
  
  # best parameter and its index
  bestind <- which.max(perf)
  
  hypergrid <- tuneres@hypergrid
  temp <- as.list(hypergrid[bestind, , drop = FALSE])
  attr(temp, "out.attrs") <- NULL
  param <- temp
  
  return(new("BestParameters", param = param, bestind = bestind, 
             perf = perf, measure = measurec))
} 
