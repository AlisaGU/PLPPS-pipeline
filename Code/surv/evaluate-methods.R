# Filename: evaluationSurv.r
# Title: Evaluation measures based on the survival results
#
# Author: Christoph Bernau 
# Email: <bernau>
# Date of creation: 17.4.2012
#
# Brief description:
#   Evaluates resampling results of different survival methods.
#
# Further comments and notes:
#   Also needed for tuning purposes.
#
# Arguments:
#   -survresult: output-object of method 'learnsurvival'
#   -measure: character indicating the evaluation metric to be used
#
# Value:
#   EvaluateOutput
###############################################################################

### signature = 'LearnOut' LearnOut is a new class containing a list of
### ModelLearned-objects and additional slots (e.g. X,y,GeneSel,nbgene) which
### could be generated at the end of the new 'classification function'
setMethod("evaluate", signature(object = "LearnOut"), function(object, 
                                                               measure = "CvPLogL", newdata = NULL, newy = NULL, 
                                                               index = NULL, imputezero = FALSE, ...) {
  
  measurec <- new(measure)
  
  ### prepare list for additional results
  survresult <- object
  additional <- list()
  ll <- list(...)
  # some administration get full X and y from LearnOut-object 
  #survresult
  
  if (is.null(newdata)) {
    Xc <- survresult@X  #getX(survresult)
    
    if (is.null(newy) != TRUE) 
      warning("Argument \"newy\" is ignored.")
    if (is.null(index) != TRUE) 
      warning("Argument \"index\" is ignored.")
    
    if (ncol(Xc) == 1 & sum(abs(Xc)) == 0) {
      Xc <- ll$X
      ll$X <- NULL
    }
    
    yc <- survresult@y  #getY(survresult)
    # get GeneSel from LearnOut
    GeneSel <- survresult@GeneSel  #getGenesel(survresult)
    # get nbgene from LearnOut
    nbgene <- survresult@nbgene  #getNbgene(survresult)
    
    # get threshold from LearnOut
    threshold <- survresult@threshold  #getNbgene(survresult)
    add <- list(...)
    foldmeasures <- list()
    
    for (i in 1:length(survresult@ModelLearnedlist)) {
      ### learnind
      learnind <- survresult@ModelLearnedlist[[i]]@learnind
      ### list for calling predict
      ll2 <- ll
      ll2$object <- survresult@ModelLearnedlist[[i]]
      ll2$type <- measurec@measuretype
      
      ### include GeneSelection results
      if (length(GeneSel@method) != 0) {
        Xi <- Xc[, GeneSel@rankings[[1]][i, 1:nbgene]]
        
        if (threshold != -1) {
          impi <- GeneSel@importance[[1]][i, ]
          
          if (GeneSel@criterion == "coefficient") 
            impi <- impi[impi > threshold] 
          else impi <- impi[impi < threshold]
          nbgene2 <- min(length(impi), nbgene)
          Xi <- Xc[, GeneSel@rankings[[1]][i, 1:nbgene2]]
          
        }
      } else {
        Xi <- Xc
      }
      
      
      
      
      # newdata
      ll2$newdata <- Xi[-learnind, ]
      
      # predict (testdata)
      predictionobject <- do.call(predict, args = ll2)
      
      # further arguments for evaluate
      newy <- yc[-learnind, ]
      ll$lptrain <- 
        survresult@ModelLearnedlist[[i]]@linear.predictor
      ll$ytrain <- yc[learnind, ]
      
      
      # evaluate
      foldmeasures[[i]] <- evaluate(predictionobject, newy, 
                                    measurec, add = ll)
    }
    out <- new("EvaluateOutput", result = foldmeasures, 
               measure = measure, 
               additional = additional)
  } else {
    if (length(survresult@ModelLearnedlist) != 1 & is.null(index) 
        == TRUE) stop("Object contains several models and 
                      no index is provided")
    
    if (is.null(index) == TRUE) 
      index <- 1
    
    if (is.null(newy)) 
      stop("Model cannot be evaluated without survival 
           response on test data. Please specify 
           \"newy\".")
    
    Xtrain <- survresult@X[survresult@
                             ModelLearnedlist[[index]]@learnind, ]
    
    ## Fill in missing variables of newdata with zero if 
    #imputezero=TRUE
    if (imputezero) {
      n.missing.vars <- sum(!colnames(Xtrain) %in% 
                              colnames(newdata))
      if (n.missing.vars > 0) {
        zeros <- matrix(0, ncol = n.missing.vars, 
                        nrow = nrow(newdata))
        colnames(zeros) <- colnames(Xtrain)[!colnames(Xtrain) 
                                            %in% colnames(newdata)]
        rownames(zeros) <- rownames(newdata)
        newdata <- cbind(newdata, zeros)
        rm(zeros)
      }
    }
    
    newdata <- newdata[, na.omit(match(colnames(Xtrain), 
                                       colnames(newdata)))]
    
    if (!identical(colnames(newdata), colnames(Xtrain))) 
      stop("newdata is missing variables that existed in the 
           training set.\nYou may either: 1) call 
           evaluate() with imputezero=TRUE, 
           or impute missing variables\nprior to 
           calling evaluate().")
    
    
    
    # prepare arguments for predict () call
    ll <- ll2 <- list(...)
    ll$measure <- measure
    ll$object <- survresult@ModelLearnedlist[[index]]
    ll$newy <- newy
    
    nbgene <- survresult@nbgene
    threshold <- survresult@threshold
    
    GeneSel <- survresult@GeneSel  #getGenesel(survresult)
    ### include GeneSelection results
    if (length(GeneSel@rankings) != 0) {
      newdata <- newdata[, GeneSel@rankings[[1]][index, 1:nbgene]]
      
      if (threshold != -1) {
        impi <- GeneSel@importance[[1]][index, ]
        
        if (GeneSel@criterion == "coefficient") 
          impi <- impi[impi > threshold] 
        else impi <- impi[impi < threshold]
        nbgene2 <- min(length(impi), nbgene)
        newdata <- newdata[, GeneSel@rankings[[1]][index, 
                                                   1:nbgene2]]
      }
      
    }
    
    
    ll$newdata <- newdata
    out <- do.call("evaluate", args = ll)
  }
  
  return(out)
})

setMethod("evaluate", signature(object = "LinearPrediction"), 
          function(object, newy, measure, add = list()) {
            
            evaluate(measure, object@lp, newy, add)
          })


### differentiation of SurvivalProbs and LinearPrediction should be maintained
setMethod("evaluate", signature(object = "SurvivalProbs"), function(object, 
                                                                    newy, measure, 
                                                                    add = list()) evaluate(measure, object@SurvivalProbs, newy,
                                                                                           add))

### signaure: ModelLearned
setMethod("evaluate", signature(object = "ModelLearned"), 
          function(object, newdata, 
                   newy, measure = "CvPLogL", add = list(), ...) {
            measurec <- new(measure)
            additional <- list()
            ### prepare arguments for predict call
            ll <- list(...)
            ll$object <- object
            ll$newdata <- newdata
            ll$type <- measurec@measuretype
            
            predictionobject <- do.call(predict, args = ll)
            # evaluate
            
            add <- c(add, list(...))
            add$ytrain <- object@y
            add$lptrain <- object@linear.predictor
            
            foldmeasures <- list(evaluate(predictionobject, newy, measurec, 
                                          add))
            new("EvaluateOutput", result = foldmeasures, measure = measure, 
                additional = additional)
          })

### signature:modelC
setMethod("evaluate", signature(object = "ModelBase"), function(object, newdata, 
                                                                newy, measure = "CvPLogL", add = list(), lptrain = NULL, 
                                                                ytrain = NULL, ...) {
  
  measurec <- new(measure)
  additional <- list()
  ### prepare arguments for predict call
  ll <- list(...)
  ll$object <- object
  ll$newdata <- newdata
  ll$type <- measurec@measuretype
  
  predictionobject <- do.call(predict, args = ll)
  # evaluate
  add <- c(add, list(...))
  add$lptrain <- lptrain
  add$ytrain <- ytrain
  
  foldmeasures <- list(evaluate(predictionobject, newy, 
                                measurec, add))
  
  new("EvaluateOutput", result = foldmeasures, measure = measure, 
      additional = additional)
})


setMethod("evaluate", signature(object = "coxph"), function(object, ...) {
  # wrap the coxph object into a ModelCustom object and then just use
  # evaluate()
  evaluate(new("ModelCustom", predfun = function(object, ...) 
    new("LinearPrediction", lp = 
          predict(object@mod, ...)), 
    mod = object, extraData = list()), ...)
})