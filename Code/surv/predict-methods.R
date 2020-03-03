setMethod("predict", signature(object = "LearnOut"), function(object, ...) {
  predictsurvhd(object, ...)
})

setMethod("predict", signature(object = "ModelLearned"), function(object, ...) {
  predictsurvhd(object, ...)
})

setMethod("predict", signature(object = "ModelBase"), function(object, ...) {
  predictsurvhd(object, ...)
})

setMethod("predict", signature(object = "NULL"), function(object, ...) {
  predictsurvhd(object, ...)
})

setGeneric("predictsurvhd", function(object, newdata, ...) 
  standardGeneric("predictsurvhd"))

setMethod("predictsurvhd", signature(object = "ANY", newdata = "ExpressionSet"), 
          function(object, newdata, ...) {
            predictsurvhd(object = object, newdata = 
                            as.data.frame(t(exprs(newdata)), 
                                          check.names = FALSE), ...)
            
          })

setMethod("predictsurvhd", signature(object = "ANY", newdata = "matrix"), 
          function(object, newdata, ...) {
            predictsurvhd(object = object, newdata = 
                            as.data.frame(newdata, check.names = FALSE), ...)
            
          })

setMethod("predictsurvhd", signature(object = "NULL", newdata = "data.frame"), 
          function(object, newdata, type, timegrid = NULL, ...) {
            pred <- rep(NA, nrow(newdata))
            names(pred) <- rownames(newdata)
            class(pred) <- "numeric"
            pred <- new("LinearPrediction", lp = pred)
            return(pred)
          })

setMethod("predictsurvhd", signature(object = "ModelLearned", 
                                     newdata = "missing"), 
          function(object, newdata, type = "lp", timegrid = NULL, 
                   gbm = TRUE, ...) {
            
            pred <- object@linear.predictor
            return(pred)
          })

setMethod("predictsurvhd", signature(object = "ModelLearned", 
                                     newdata = "data.frame"), 
          function(object, newdata, type = "lp", timegrid = NULL, 
                   gbm = TRUE, ...) {
            
            ll <- list(...)
            if (type == "lp") {
              if (!gbm) 
                warning("Argument gbm=FALSE ignored since prediction 
                        type is \"lp\".")
              pred <- predict(object@model, type = type, 
                              newdata = newdata, ll)
            } else if (type == "SurvivalProbs") {
              if (gbm) {
                survbase <- gbmbasehaz(object, timegrid = timegrid)
                newlp <- predict(object@model, newdata, type = "lp", ll)@lp
                ## individual survival probabilities
                curves <- matrix(survbase, ncol = length(timegrid), 
                                 nrow = length(newlp), 
                                 byrow = TRUE)^(exp(newlp))
                pred <- new("breslow", curves = curves, time = timegrid)
                pred <- new("SurvivalProbs", SurvivalProbs = pred)
              } else {
                pred <- predict(object@model, newdata = newdata, type = 
                                  type, timegrid = timegrid, ll)
              }
            }
            return(pred)
          })



setMethod("predictsurvhd", signature(object = "LearnOut", newdata = "missing"), 
          function(object, newdata, type = "lp", timegrid = NULL, gbm = TRUE, 
                   indices = NULL, addtrain = FALSE, 
                   voting_scheme = NULL, ...) {
            ll <- list(...)
            if (ncol(object@X) == 1 & sum(abs(object@X)) == 0) 
              stop("No new data specified and no data stored in slot \"X\"
                   of the LearnOut-object.")
            
            if (type != "cp") {
              
              ## prepare argument
              ll$gbm <- gbm
              orgX <- object@X
              ll$type <- type
              ll$timegrid <- timegrid
              
              pred <- list()
              
              if (is.null(indices)) 
                indices <- 1:length(object@ModelLearnedlist)
              for (j in indices) {
                ll$object <- object@ModelLearnedlist[[j]]
                
                if (!addtrain) 
                  ll$newdata <- orgX[-object@ModelLearnedlist[[j]]
                                     @learnind, ] else ll$newdata <- orgX
                                     
                                     predj <- do.call("predict", args = ll)
                                     
                                     if (type == "lp") {
                                       names(predj@lp) <- rownames(ll$newdata)
                                     } else if (type == "SurvivalProbs") {
                                       rownames(predj@SurvivalProbs@curves) <- 
                                         rownames(ll$newdata)
                                     }
                                     
                                     pred[[j]] <- predj
              }
            } else {
              orgX <- object@X
              
              if (addtrain) 
                warning("Parameter \"addtrain\" ignored.")
              
              if (is.null(ll$method)) 
                ll$method <- "median"
              
              # first calculate the cutpoints for each fold in the 
              #learnind samples only
              cutpoints = lapply(object@ModelLearnedlist, function(xx)
                plotKMStratifyBy(method = ll$method, 
                                 y = xx@y, linearriskscore = 
                                   xx@linear.predictor@
                                   lp, labels = ll$labels, 
                                 plot = FALSE))
              
              # now do the patient stratification in the validation 
              #folds /validation set
              res <- as.factor(do.call(c, lapply(1:length(cutpoints), 
                                                 function(i) {
                                                   newdat <- orgX[-object@ModelLearnedlist[[i]]@learnind, ]
                                                   
                                                   # get the risk scores
                                                   lrs <- predict(object@ModelLearnedlist[[i]], 
                                                                  newdata = newdat)@lp
                                                   # do the stratification
                                                   strata <- .kmStratify(NULL, linearriskscore = lrs, 
                                                                         cutpoints = cutpoints[[i]]$cutpoints, 
                                                                         plot = FALSE, labels=ll$labels)$strata
                                                   # this ugliness was the shortest way I figured 
                                                   #out to combine factors with the
                                                   # same levels in R. Gosh.
                                                   strata <- as.character(strata)
                                                   names(strata) <- names(lrs)
                                                   strata
                                                 })))
              ###
              allobs <- rownames(orgX)
              
              clmat <- matrix(ncol = max(table(names(res))),
                              nrow = length(allobs), NA, 
                              dimnames = list(allobs, 1:max(table(names(res)))))
              for (cc in 1:length(allobs)) {
                rescc <- res[which(names(res) == allobs[cc])]
                clmat[cc, ] <- c(as.character(rescc),
                                 rep(NA, times = ncol(clmat) - length(rescc)))
              }
              
              if (!is.null(voting_scheme)) {
                if (class(voting_scheme) != "function") {
                  
                  if (voting_scheme == "first") 
                    pred <- apply(clmat, 1, FUN = function(x) {
                      names(table(x))[which.max(table(x))]
                    }) else if (voting_scheme == "last") 
                      pred <- apply(clmat, 1, FUN = function(x) {
                        names(table(x))[max(which(table(x) 
                                                  == max(table(x))))]
                      }) else if (voting_scheme == "random") 
                        pred <- apply(clmat, 1, FUN = function(x) {
                          names(table(x))[sample(which(table(x) 
                                                       == max(table(x))), 1)]
                        })
                      
                } else pred <- apply(clmat, 1, FUN = voting_scheme)
              } else pred <- clmat
            }
            return(pred)
          })


setMethod("predictsurvhd", signature(object = "LearnOut", 
                                     newdata = "data.frame"), 
          function(object, newdata, type = "lp", timegrid = NULL, 
                   gbm = TRUE, indices = NULL, addtrain = FALSE, 
                   voting_scheme = NULL, ...) {
            ll <- list(...)
            
            if (type != "cp") {
              ## prepare argument
              ll$gbm <- gbm
              ll$type <- type
              ll$timegrid <- timegrid
              
              pred <- list()
              
              if (is.null(indices)) 
                indices <- 1:length(object@ModelLearnedlist)
              
              for (j in indices) {
                ll$object <- object@ModelLearnedlist[[j]]
                
                if (addtrain) {
                  if (ncol(object@X) == 1 & sum(abs(object@X)) == 0) 
                    stop("No data stored in slot \"X\" of the 
                         LearnOut-object, cannot add training data!")
                  ll$newdata <- rbind(object@X, newdata)
                } else (!addtrain)
                ll$newdata <- newdata
                
                predj <- do.call("predict", args = ll)
                
                if (type == "lp") {
                  names(predj@lp) <- rownames(ll$newdata)
                } else if (type == "SurvivalProbs") {
                  rownames(predj@SurvivalProbs@curves) <- 
                    rownames(ll$newdata)
                }
                
                pred[[j]] <- predj
              }
              } else {
                if (addtrain) 
                  warning("Parameter \"addtrain\" ignored.")
                
                if (is.null(ll$method)) 
                  ll$method <- "median"
                
                # first calculate the cutpoints for each fold in the learnind 
                #samples only
                cutpoints = lapply(object@ModelLearnedlist, function(xx)
                  plotKMStratifyBy(method = ll$method, 
                                   y = xx@y, linearriskscore = 
                                     xx@linear.predictor@lp, plot = FALSE))
                
                # now do the patient stratification in the validation
                # folds /validation set
                res <- as.factor(do.call(c, lapply(1:length(cutpoints), 
                                                   function(i) {
                                                     # get the risk scores
                                                     lrs <- predict(object@ModelLearnedlist[[i]], 
                                                                    newdata = newdata)@lp
                                                     # do the stratification
                                                     strata <- .kmStratify(NULL, linearriskscore = 
                                                                             lrs, cutpoints = cutpoints[[i]]$cutpoints, 
                                                                           plot = FALSE, labels= ll$labels)$strata
                                                     # this ugliness was the shortest way 
                                                     #I figured out to combine factors with the
                                                     # same levels in R. Gosh.
                                                     strata <- as.character(strata)
                                                     names(strata) <- names(lrs)
                                                     strata
                                                   })))
                allobs <- rownames(newdata)
                
                clmat <- matrix(ncol = max(table(names(res))), 
                                nrow = length(allobs), 
                                NA, dimnames = list(allobs, 1:max(table(names(res)))))
                for (cc in 1:length(allobs)) {
                  rescc <- res[which(names(res) == allobs[cc])]
                  clmat[cc, ] <- c(as.character(rescc), rep(NA, 
                                                            times = ncol(clmat) - 
                                                              length(rescc)))
                }
                
                if (!is.null(voting_scheme)) {
                  if (class(voting_scheme) != "function") {
                    if (voting_scheme == "first") 
                      pred <- apply(clmat, 1, FUN = function(x) {
                        names(table(x))[which.max(table(x))]
                      }) else if (voting_scheme == "last") 
                        pred <- apply(clmat, 1, FUN = function(x) {
                          names(table(x))[max(which(table(x)
                                                    == max(table(x))))]
                        }) else if (voting_scheme == "random") 
                          pred <- apply(clmat, 1, FUN = function(x) {
                            names(table(x))[sample(which(table(x)
                                                         == max(table(x))), 1)]
                          })
                        
                  } else pred <- apply(clmat, 1, FUN = voting_scheme)
                } else pred <- clmat
              }
            
            return(pred)
          })