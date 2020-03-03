# Filename:customSurv.r
# Title: integrates user defined surviva methods into learnSurvival.
#
# Author: Christoph Bernau, adapted from M. Slawski and A-L Boulesteix
# Email: <bernau@ibe.med.uni-muenchen.de>
# Date of creation: April 5, 2012
#
# Brief description:
#   Returns an object of class ModelLearned.
#
# Arguments:
#   -X: matrix of variables (rows are observations,columns are variables)
#   -y: survival response of class Surv 
#   -learnind: vector indicating which observations are to be used for
#   training 
#              the survival model
#    -customSurvModel: function accepting X,y, learnind and ... and returns an
#    object of class ModelCustom
#
# Value:
#   ModelLearned
#
###############################################################################



setGeneric("customSurv", function(X, y, ...) standardGeneric("customSurv"))

setMethod("customSurv", signature(X = "data.frame", y = "Surv"), function(X, y,
                                                                          learnind, customSurvModel, essential = FALSE, ...) {
  
  if (class(customSurvModel) != "function") 
    stop(" Argument \"customSurvModel\" must be of class 
         \"function\".")
  
  nrx <- nrow(X)
  yAll <- y
  if (nrx != nrow(y)) 
    stop("Number of rows of 'X' must agree with length of y \n")
  if (missing(learnind)) 
    learnind <- 1:nrx
  
  if (length(learnind) > nrx) 
    stop("length of 'learnind' must be smaller than the number of 
         observations. \n")
  
  mode <- "learnsurvival"
  
  if (!identical(all.equal(colnames(X), make.names(colnames(X), 
                                                   unique = TRUE)), TRUE)) {
    colnames(X) <- make.names(colnames(X), unique = TRUE)
    # warning('converting colnames(X) to valid, unique names')
  }
  
  Ylearn <- y[learnind]
  Xlearn <- data.frame(X[learnind, ])
  # dotsCall <- list(...) ll <- eval(dotsCall)
  ll <- list(...)
  
  ll$Ylearn <- Ylearn
  ll$learnind <- learnind
  ll$Xlearn <- Xlearn
  
  output.custom <- do.call("customSurvModel", args = ll)
  
  Xtest <- data.frame(X[-learnind, ])
  
  if (nrow(Xtest) == 0) {
    Xtest <- Xlearn
    y <- Ylearn
  } else y <- y[-learnind]
  
  pred.learn <- try(predict(output.custom, newdata = Xlearn,
                            type = "lp"), silent = TRUE)
  
  if (class(pred.learn) == "try-error") 
    pred.learn <- new("LinearPrediction")
  
  if (essential == TRUE) {
    2 + 2
    # prune modd
  }
  
  return(new("ModelLearned", y = Ylearn, linear.predictor = 
               pred.learn, learnind = learnind, 
             method = "customSurv", model = output.custom))
})

setMethod("customSurv", signature(X = "matrix", y = "Surv"), 
          function(X, y, ...) {customSurv(X = as.data.frame(X, 
                                                            check.names = FALSE), y = y, ...)})

setMethod("customSurv", signature(X = "ExpressionSet", y = "Surv"), 
          function(X, y, ...) { customSurv(X = as.data.frame(t(exprs(X))
                                                             , check.names = FALSE), y = y, ...)})

setMethod("customSurv", signature(X = "ExpressionSet", y = "character"), 
          function(X, y, ...) {
            Xdat <- as.data.frame(t(exprs(X)), check.names = FALSE)
            customSurv(X = Xdat, y = .fetchyFromEset(X, y), ...)
          })

setMethod("predictsurvhd", signature(object = "ModelCustom", 
                                     newdata = "data.frame"), function(object, newdata, ...) {
                                       funargs <- list(...)
                                       funargs$object <- object
                                       funargs$newdata <- newdata
                                       pred <- do.call(object@predfun, args = funargs)
                                       return(pred)
                                     }) 
