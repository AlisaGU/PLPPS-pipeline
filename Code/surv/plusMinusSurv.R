# Filename: plusminusCMA.r
# Title: One of many classifiers.
#
# Author: Levi Waldron, adapted from M. Slawski and A-L Boulesteix
# Email: <lwaldron.research@gmail.com>
# Date of creation: Mar 21, 2012
#
# Brief description:
#   Returns an object of class ModelLearned.
#
# Further comments and notes:
#   The plusminus learner assigns coefficients of +/-1 to each feature,
#   optionally scaled by the standard deviation of the feature for
#   equal contributions to the continuous risk score.  Signs are
#   determined from the sign of the univariate Cox regression for each
#   feature.
#
# Arguments:
#   -X: matrix of variables (rows are observations,columns are variables)
#   -y: survival response of class Surv 
#   -learnind: vector indicating which observations are to be used for 
#                training the survival model
#   -modeltype: "plusminus" for +1/-1 coefficients with optional
#               scaling, compoundcovariate for coefficients equal to the Cox
#               coefficient, tscore for the TCGA-type model (score is statistic 
#                for t-test between good and bad-prognosis genes).
#
# Value:
#   ModelLearned
#
###############################################################################

setGeneric("plusMinusSurv", function(X, y, ...) 
  standardGeneric("plusMinusSurv"))

setMethod("plusMinusSurv", signature(X = "data.frame", y = "Surv"), 
          function(X, y, learnind, essential = FALSE,
                   modeltype = "plusminus", ...) {
            
            nrx <- nrow(X)
            yAll <- y
            if (nrx != nrow(y)) 
              stop("Number of rows of 'X' must agree with length of y \n")
            if (missing(learnind)) 
              learnind <- 1:nrx
            
            if (length(learnind) > nrx) 
              stop("length of 'learnind' must be smaller than the number of
                   observations. \n")
            
            Xlearn <- data.frame(X[learnind, ], check.names = FALSE)
            Ylearn <- y[learnind]
            
            ll <- list(...) 
            #test this alternative to eval(substitute(list(...)))
            
            ll$X <- Xlearn
            ll$y <- Ylearn
            ll$modeltype <- modeltype
            
            good.arg.names <- names(formals(plusMinus))
            output <- do.call("plusMinus", args = ll[names(ll)
                                                     %in% good.arg.names])
            
            ## Predictions on learning set can be used later for 
            ##calculating baseline hazards
            linear.predictor <- predict(output, newdata = Xlearn, type = "lp")
            
            if (essential == TRUE) {
              2 + 2
              # prune output here
            }
            
            ## ModelLearned should contain only method, list of the core object
            ##, learnind.
            return( new("ModelLearned", y = Ylearn, linear.predictor = 
                          linear.predictor, learnind = learnind,
                        method = "plusMinusSurv", model = output) )
          })

setMethod("plusMinusSurv", signature(X = "matrix", y = "Surv"),
          function(X, y, ...) {
            plusMinusSurv(X = as.data.frame(X, check.names = FALSE), y = y, ...)
          })

setMethod("plusMinusSurv", signature(X = "ExpressionSet", y = "Surv"), 
          function(X, y, ...) {
            plusMinusSurv(X = as.data.frame(t(exprs(X)), 
                                            check.names = FALSE), y = y, ...)
          })

setMethod("plusMinusSurv", signature(X = "ExpressionSet", y = "character"),
          function(X, y, ...) {
            Xdat <- as.data.frame(t(exprs(X)), check.names = FALSE)
            plusMinusSurv(X = Xdat, y = .fetchyFromEset(X,y), ...)
          })

setMethod("predictsurvhd", signature(object = "ModelLinear",
                                     newdata = "data.frame"), 
          function(object, newdata, type="lp", timegrid = NULL, ...) {
            ll2 <- list(...)
            dropmissing <- ll2$ll$dropmissing
            if (is.null(dropmissing))
              dropmissing <- TRUE
            
            if (!length(object@modeltype)) object@modeltype = 
              "compoundcovariate"
            
            if (type == "lp") {
              newdata <- as.matrix(newdata)
              cc <- object@coefficients
              ## does the model have an intercept?
              if (grepl("(Intercept)", names(cc)[1], fixed = TRUE)) {
                intercept <- cc[1]
                cc <- cc[-1]
              } else {
                intercept <- 0
              }
              if (dropmissing) {
                cc <- cc[names(cc) %in% colnames(newdata)]
              }
              rowN<-rownames(newdata)
              colN<-names(cc)
              #newdata <- newdata[, colnames(newdata) %in% names(cc)]
              newdata<-matrix(newdata[, colnames(newdata) %in% names(cc)],
                              ncol = length(which(colnames(newdata) %in% names(cc))))
              rownames(newdata)<-rowN
              colnames(newdata)<-colN
              cc <- cc[match(colnames(newdata), names(cc))]
              if (!identical(names(cc), colnames(newdata))) {
                ## can happen if none of the coefficients are found 
                ##in newdata
                warning("names of coefficients do not match colnames
                        (newdata).  Make sure columns of newdata 
                        correspond to variables.")
                pred <- rep(NA, nrow(newdata))
                names(pred) <- rownames(newdata)
              } else {
                if (object@modeltype == "tscore") {
                  require(genefilter)
                  used.features <- abs(cc) > 0
                  cc <- cc[used.features]
                  newdata <- newdata[, used.features]
                  if (length(unique(sign(cc))) == 1) 
                    warning("Must have at least one good-prognosis and 
                            one bad-prognosis feature to make 
                            predictions 
                            with the tscore method.  
                            Returning all NA predictions.")
                  pred <- try(-rowttests(newdata, fac = factor(sign(cc)), 
                                         tstatOnly = TRUE)[, 1])
                } else if (object@modeltype %in% c("voting", 
                                                   "positiveriskvoting", 
                                                   "negativeriskvoting")) {
                  votingthresholds <- object@votingthresholds
                  if (length(votingthresholds) == 1) 
                    votingthresholds <- structure(rep(votingthresholds,
                                                      length(cc)), .Names = names(cc))
                  votingthresholds <- votingthresholds[match(names(cc), 
                                                             names(votingthresholds))]
                  if (!identical(all.equal(names(votingthresholds), 
                                           names(cc)), TRUE)) 
                    stop("length(model@votingthresholds) must be 1 or 
                         length(cc).  If length(cc), their 
                         names must match.")
                  newdata.dichot <- sweep(newdata, 2, votingthresholds)
                  ## +1 if greather than threshold, -1 otherwise
                  newdata.dichot <- sign(newdata.dichot)
                  ## multiply by coefficints.  Value above threshold, 
                  ##positive coef -> positive
                  ## contribution to risk, etc.
                  newdata.dichot <- sweep(newdata.dichot, 2, cc, "*")
                  # reset any negative/positive votes to zero
                  if (object@modeltype == "positiveriskvoting") 
                    newdata.dichot[newdata.dichot < 0] <- 0 else if 
                  (object@modeltype == "negativeriskvoting") 
                      newdata.dichot[newdata.dichot > 0] <- 0
                  ## predictions are just the sum of the votes:
                  pred <- rowSums(newdata.dichot)
                } else {
                  pred <- try((as.matrix(newdata) %*% cc)[, 1]
                              + intercept)
                }
                if (class(pred) == "try-error") {
                  pred <- rep(NA, nrow(newdata))
                  names(pred) <- rownames(newdata)
                } else {
                  pred <- structure(pred, .Names = rownames(pred))
                }
              }
              if (length(object@cutpoint) > 0) {
                pred <- cut(pred, breaks = c(-Inf, object@cutpoint, Inf), 
                            ordered_result = TRUE, 
                            include.lowest = TRUE)
                pred <- as.numeric(pred)
              }
              pred <- new("LinearPrediction", lp = pred)
            } else if (type == "SurvivalProbs") {
              stop('Survival method does not provide an own implementation 
                   survival probability predicition. \n You could 
                   try "gbm=TRUE".')
              # warning("Coxbased survival-probability estimation
              #(gbm set to TRUE).")
              # predict(object, newdata, type = "SurvivalProbs",
              #timegrid, gbm = TRUE)
            }
            return(pred)
            })

