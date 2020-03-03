# Filename: coxBoostSurv.r
# Title: One of many survival methods.
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
#   -learnind: vector indicating which observations are to be used for training 
#              the survival model
#   -unpenalizedvars: matrix of unpenalized variables
#
# Value:
#   ModelLearned
#
###############################################################################



setGeneric("coxBoostSurv", function(X, y, ...) standardGeneric("coxBoostSurv"))

setMethod("coxBoostSurv", signature(X = "data.frame", y = "Surv"), 
          function(X, y, learnind, unpenalizedvars = NULL, 
                   essential = FALSE,...) {
            require(CoxBoost)
            
            if (!is.null(unpenalizedvars)) 
              unpeninds <- 1:ncol(unpenalizedvars)
            
            nrx <- nrow(X)
            yAll <- y
            if (nrx != nrow(y)) 
              stop("Number of rows of 'X' must agree with length of y \n")
            if (missing(learnind)) 
              learnind <- 1:nrx
            
            if (length(learnind) > nrx) 
              stop("length of 'learnind' must be smaller than the number 
                   of observations. \n")
            
            mode <- "learnsurvival"
            
            if (!identical(all.equal(colnames(X), make.names(colnames(X), 
                                                             unique = TRUE)), TRUE)) {
              colnames(X) <- make.names(colnames(X), unique = TRUE)
              # warning('converting colnames(X) to valid, unique names')
            }
            
            Ylearn <- y[learnind]
            
            if (is.null(unpenalizedvars) == TRUE) 
              Xlearn <- data.frame(X[learnind, ])
            else Xlearn <- data.frame(unpenalizedvars[learnind, ], 
                                      X[learnind, ])
            
            ll <- list(...)
            
            ### check that user does not try to specify unpen.index
            if (is.element("unpen.index", names(ll))) 
              stop("Please use the argument \"unpenalizedvars\" for 
                   specifying mandatory variables. ")
            
            ll$time <- Ylearn[, 1]
            ll$status <- Ylearn[, 2]
            ll$trace = FALSE
            ll$standardize = TRUE
            if (!is.null(unpenalizedvars)) 
              ll$unpen.index <- unpeninds
            
            ll$x <- as.matrix(Xlearn)
            
            output.coxboost <- do.call("CoxBoost", args = ll)
            
            if (is.null(unpenalizedvars)) Xtest <- data.frame(X[-learnind, ])
            else Xtest <- data.frame(unpenalizedvars[-learnind, ],
                                     X[-learnind, ])
            
            if (nrow(Xtest) == 0) {
              Xtest <- Xlearn
              y <- Ylearn
            } else y <- y[-learnind]
            #####
            
            output.coxboost <- new("CoxBoostSurv", mod = output.coxboost)
            
            
            ## Linear risk score
            
            pred.learn <- predict(output.coxboost, newdata = Xlearn, 
                                  type = "lp")
            
            if (essential == TRUE) {
              2 + 2
              # prune modd
            }
            
            return( new("ModelLearned", y = Ylearn, 
                        linear.predictor = pred.learn, 
                        learnind = learnind, method = "coxBoostSurv", 
                        model = output.coxboost) )
          })

setMethod("coxBoostSurv", signature(X = "matrix", y = "Surv"), function(X, y,
                                                                        ... ) {
  coxBoostSurv(X = as.data.frame(X, check.names = FALSE), 
               y = y, ... )
})

setMethod("coxBoostSurv", signature(X = "ExpressionSet", y = "Surv"), 
          function(X, y,... ) {
            coxBoostSurv(X = as.data.frame(t(exprs(X)), check.names = FALSE), 
                         y = y, ... )
          })

setMethod("coxBoostSurv", signature(X = "ExpressionSet", y = "character"), 
          function(X, y, ... ) {
            Xdat <- as.data.frame(t(exprs(X)), check.names = FALSE)
            coxBoostSurv(X = Xdat, y = .fetchyFromEset(X, y), ... )
          })

setMethod("predictsurvhd", signature(object = "CoxBoostSurv", 
                                     newdata = "data.frame"), 
          function(object, newdata, type, timegrid = NULL, ...) {
            require(CoxBoost)
            if (type == "lp") {
              survobj <- object@mod
              pred <- predict(object = survobj, type = "lp", 
                              newdata = newdata)
              pred <- structure(pred[1, ], .Names = colnames(pred))
              pred <- new("LinearPrediction", lp = pred)
            } else if (type == "SurvivalProbs") {
              survobj <- object@mod
              if (is.null(timegrid)) {
                stop("No timegrid specified.")
              }
              curves <- predict(object = survobj, newdata = newdata, 
                                times = timegrid, type = "risk")
              pred <- new("breslow", curves = curves, time = timegrid)
              pred <- new("SurvivalProbs", SurvivalProbs = pred)
            }
            
            else stop('Invalid "type" argument.')
            
            return(pred)
          })
