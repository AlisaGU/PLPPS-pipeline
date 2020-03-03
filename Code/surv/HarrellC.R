# Filename: harrellC.r
# Title: Determine harrell's C on a test set.
#
# Author: Christoph Bernau
# Email: <bernau@ibe.med.uni-muenchen.de>
# Date of creation: April 27, 2012
#
# Brief description:
#   Returns value of harrell's C.
#
# Arguments:
#   -pred.test: the linear predictor fr the new test observations
#   -newy: the survival response for the test cases (class Surv)
#
# Value:
#   value of Harrell's C-statistic
#
###############################################################################




setMethod("evaluate", signature(object = "HarrellC"), function(object, predtest, 
                                                               newy, add = list(),...) {
  
  if (!class(predtest) %in% c("numeric", "integer")) 
    stop("harrellC: predtest must be numeric or integer")
  if (!class(newy) %in% "Surv") 
    stop("harrellC: newy must be of class Surv")
  ## done argument checking
  harrell.C <- try(rcorr.cens(x = -predtest, S = newy)["C Index"], 
                   silent = FALSE)
  if (class(harrell.C) == "try-error") 
    harrell.C <- NA
  return(harrell.C)
})

setMethod("aggregateMeasure", signature(object = "HarrellC"), 
          function(object, foldmeasures) {
            return(mean(unlist(foldmeasures)))
          }) 