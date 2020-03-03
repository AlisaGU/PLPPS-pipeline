# Filename: UnoC.r
# Title: Determine Uno's C on a test set.
#
# Author: Christoph Bernau
# Email: <bernau@ibe.med.uni-muenchen.de>
# Date of creation: April 27, 2012
#
# Brief description:
#   Returns value of Uno's C.
#
# Further comments and notes:
#   
# Arguments:
#   -predtest: linear predictor for test cases
#   -newy: survival response for test cases (class Surv)
#   -tau: truncation time
#
# Value:
#   value of Uno's C-statistic
#
######################################################################


setMethod("evaluate", signature = (object = "UnoC"), function(object, predtest,
                                                              newy, add = list(),...) {
  
  tau <- add$tau
  ## make sure to set seed here, or Inf.Cval will set
  ##and synchronize it!
  if (is.null(add$seed)) add$seed <- sample(1e+09, 1) 
  if (!is.null(add$itr)) 
    warning(paste("Confidence interval estimation currently not",
                  "supported. Ignoring itr argument."))
  add$itr  <- 1
  # predtest,newy,tau
  if (!class(predtest) %in% c("numeric", "integer")) 
    stop("unoC: predtest must be numeric or integer")
  if (!class(newy) %in% "Surv") 
    stop("unoC: newy must be of class Surv")
  if (is.null(tau)) 
    stop("Please specify tau when calculating Uno's C-statistic 
         (see ?survC1::Inf.Cval)")
  if (!class(tau) %in% c("numeric", "integer")) 
    stop("unoC: tau must be numeric or integer")
  if (tau < min(newy[, 1])) 
    stop("Tau must be greater than the smallest event time.")
  ## setup
  cens <- try(kmcens(newy[, 1], newy[, 2], tau), silent=FALSE)
  if(class(cens) == "try-error"){
    uno.C <- NA
  }else{
    ##kmcens from survC1 package
    GXi <- cens$surv[match(newy[, 1], cens$distinct, nomatch = 1)]
    Wi <- 1/GXi/GXi * newy[, 2] * as.numeric(newy[, 1] < tau)
    uno.C <- conc(newy[, 1], newy[, 2], Wi, predtest)
  }
  return(uno.C)
})

setMethod("aggregateMeasure", signature(object = "UnoC"), 
          function(object, foldmeasures) { return(mean(unlist(foldmeasures)))
          }) 