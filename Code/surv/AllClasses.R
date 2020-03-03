###############################################################################
###############################################################################
###old classes
setOldClass("Surv")


###############################################################################
###############################################################################
###basic model classes
###superclass
setClass(Class = "ModelBase", representation(mod = "ANY"))

##accessor(s)
#generic(s)
setGeneric("mod", function(object) standardGeneric("mod"))


#method(s)
setMethod("mod", signature = "ModelBase", definition = 
            function(object) object@mod)


##replacement function(s)
#generic(s)
setGeneric("mod<-", function(object, value) standardGeneric("mod<-"))


#method(s)
setMethod("mod<-", signature = c("ModelBase", "ANY"), 
          definition = function(object, value){
            ##checks??  
            object@mod <- value
            object
          })




### embedding-class for coxboostsurv (needed for predict) (no additional slots)
setClass("CoxBoostSurv", contains = "ModelBase")


###############################################################################
####extended model classes
### embedding-class for custom models (needed for predict)
setClass(Class = "ModelCustom", representation(predfun = "function", 
                                               extraData = "list"), contains = "ModelBase")

##accessor(s)
#generic(s)
setGeneric("predfun", function(object) standardGeneric("predfun"))
setGeneric("extraData", function(object) standardGeneric("extraData"))

#method(s)
setMethod("predfun", signature = "ModelCustom", definition = 
            function(object) object@predfun)

setMethod("extraData", signature = "ModelCustom", definition = 
            function(object) object@extraData)


##replacement function(s)
#generic(s)
setGeneric("predfun<-", function(object, value) standardGeneric("predfun<-"))
setGeneric("extraData<-", function(object, value) 
  standardGeneric("extraData<-"))

#method(s)
setMethod("predfun<-", signature = c("ModelCustom", "function"), 
          definition = function(object, value){
            ##checks??  
            if(any( c("object","newdata","type") %in% names(formals(pf)))
               ==FALSE)
              stop('Argument names "object","newdata" and "type" 
                   are obligatory for function "prefun".')
            object@predfun <- value
            object
          })

setMethod("extraData<-", signature = c("ModelCustom", "list"), 
          definition = function(object, value){
            ##checks??  
            object@extraData <- value
            object
          })



###model linear
setClass(Class = "ModelLinear", representation(coefficients = "numeric", 
                                               cutpoint = "numeric", modeltype = "character", 
                                               votingthresholds = "numeric"), contains = "ModelBase")


##accessor(s)
#generic(s)
#setGeneric("coefficients", function(object) standardGeneric("coefficients"))


#method(s)
setMethod("coefficients", signature = "ModelLinear", definition = 
            function(object) object@coefficients)


##replacement function(s)
#generic(s)
setGeneric("coefficients<-", function(object, value) 
  standardGeneric("coefficients<-"))


#method(s)
setMethod("coefficients<-", signature = c("ModelLinear", "numeric"), 
          definition = function(object, value){
            ##checks??  
            object@coefficients <- value
            object
          })

###Constructor
ModelLinear<-function(coefficients, cutpoint, modeltype, votingthresholds){
  ###insert check of inputs
  if(!is.numeric(coefficients))
    stop("coefficients should be a numeric vector")
  if(!( length(coefficients) > 0))
    stop("coefficients should have length > 0")
  if(!is.null(cutpoint))
    if(!is.numeric(cutpoint) | !length(cutpoint)==1)
      stop("cutpoint should be a numeric vector of length 1")
  supported.modeltypes <- c("plusminus", "compoundcovariate", "tscore", 
                            "voting", "positiveriskvoting", "negativeriskvoting")
  if(!is.character(modeltype) |
     !length(modeltype)==1 |
     !modeltype[1] %in% supported.modeltypes){
    stop(paste("modeltype must be a character vector of length 1,
               and be one of the supported modeltypes: ",
               paste(supported.modeltypes, collapse=", "), sep=""))
  }
  if(!is.null(votingthresholds)){
    if(!identical(length(votingthresholds), 1)){
      if(!identical(all.equal(names(votingthresholds), names(cc)), TRUE)) 
        stop("length(votingthresholds) must be 1 or
             length(coefficients).  If length(coefficients), their names must match.")
    }}
  obj<-new("ModelLinear", coefficients=coefficients, cutpoint=cutpoint, modeltype=modeltype,
           votingthresholds=votingthresholds)
  return(obj)
    }




###############################################################################
###############################################################################
####classes necessary for evaluation
###linear prediction
setClass(Class = "LinearPrediction", representation(lp = "numeric"))

##accessor(s)
#generic(s)
setGeneric("lp", function(object) standardGeneric("lp"))


#method(s)
setMethod("lp", signature = "LinearPrediction", definition = 
            function(object) object@lp)


##replacement function(s)
#generic(s)
setGeneric("lp<-", function(object, value) standardGeneric("lp<-"))


#method(s)
setMethod("lp<-", signature = c("LinearPrediction", "numeric"), 
          definition = function(object, value){
            ##checks??  
            object@lp <- value
            object
          })


###survival probs
setClass(Class = "SurvivalProbs", representation(SurvivalProbs = "breslow"))

##accessor(s)
#generic(s)
setGeneric("SurvivalProbs", function(object) standardGeneric("SurvivalProbs"))


#method(s)
setMethod("SurvivalProbs", signature = "SurvivalProbs", definition = 
            function(object) object@SurvivalProbs)


##replacement function(s)
#generic(s)
setGeneric("SurvivalProbs<-", function(object, value) 
  standardGeneric("SurvivalProbs<-"))


#method(s)
setMethod("SurvivalProbs<-", signature = c("SurvivalProbs", "breslow"), 
          definition = function(object, value){
            ##checks?? (s. class breslow)
            object@SurvivalProbs <- value
            object
          })




###############################################################################
###############################################################################

setClass(Class = "ModelLearned", representation(y = "Surv", 
                                                linear.predictor = "LinearPrediction", learnind = "numeric", 
                                                method = "character", model = "ModelBase"))

##accessor(s)
#generic(s)
setGeneric("y", function(object) standardGeneric("y"))
setGeneric("linear.predictor", function(object) 
  standardGeneric("linear.predictor"))
setGeneric("learnind", function(object) standardGeneric("learnind"))
setGeneric("method", function(object) standardGeneric("method"))
setGeneric("model", function(object) standardGeneric("model"))

#method(s)
setMethod("y", signature = "ModelLearned", definition = 
            function(object) object@y)
setMethod("linear.predictor", signature = "ModelLearned", definition = 
            function(object) object@linear.predictor)
setMethod("learnind", signature = "ModelLearned", definition = 
            function(object) object@learnind)
setMethod("method", signature = "ModelLearned", definition = 
            function(object) object@method)
setMethod("model", signature = "ModelLearned", definition = 
            function(object) object@model)

##replacement function(s)
#generic(s)
setGeneric("y<-", function(object, value) standardGeneric("y<-"))
setGeneric("linear.predictor<-", function(object, value) 
  standardGeneric("linear.predictor<-"))
setGeneric("learnind<-", function(object, value) standardGeneric("learnind<-"))
setGeneric("method<-", function(object, value) standardGeneric("method<-"))
setGeneric("model<-", function(object, value) standardGeneric("model<-"))


#method(s)
setMethod("y<-", signature = c("ModelLearned", "Surv"), 
          definition = function(object, value){
            ##checks??  
            object@y <- value
            object
          })

setMethod("linear.predictor<-", signature =
            c("ModelLearned", "LinearPrediction"), 
          definition = function(object, value){
            ##checks-->s. class LinearPrediction  
            object@linear.predictor <- value
            object
          })

setMethod("learnind<-", signature = c("ModelLearned", "numeric"), 
          definition = function(object, value){
            ##checks 
            if(any(value<0))
              stop('Negative value in "learnind".')
            object@learnind <- value
            object
          })

setMethod("method<-", signature = c("ModelLearned", "character"), 
          definition = function(object, value){
            if(value %in%c('plusMinusSurv','coxBoostSurv','customSurv')==FALSE)
              stop('Invalid "method" name.')
            object@method <- value
            object
          })

setMethod("model<-", signature = c("ModelLearned", "ModelBase"), 
          definition = function(object, value){
            ##checks->s. class ModelBase
            object@model <- value
            object
          })

###############################################################################
###############################################################################
###learningsets class
setClass(Class = "LearningSets", representation(learnmatrix = "matrix", 
                                                method = "character", ntrain = "numeric", iter = "numeric"))

##accessor(s)
#generic(s)
setGeneric("learnmatrix", function(object) standardGeneric("learnmatrix"))
setGeneric("method", function(object) standardGeneric("method"))
setGeneric("ntrain", function(object) standardGeneric("ntrain"))
setGeneric("iter", function(object) standardGeneric("iter"))


#method(s)
setMethod("learnmatrix", signature = "LearningSets", definition = 
            function(object) object@learnmatrix)
setMethod("method", signature = "LearningSets", definition = 
            function(object) object@method)
setMethod("ntrain", signature = "LearningSets", definition = 
            function(object) object@ntrain)
setMethod("iter", signature = "LearningSets", definition = 
            function(object) object@iter)

##replacement function(s)
#generic(s)
setGeneric("learnmatrix<-", function(object, value) 
  standardGeneric("learnmatrix<-"))
setGeneric("method<-", function(object, value) 
  standardGeneric("method<-"))
setGeneric("ntrain<-", function(object, value) 
  standardGeneric("ntrain<-"))
setGeneric("iter<-", function(object, value) 
  standardGeneric("iter<-"))

#method(s)
setMethod("learnmatrix<-", signature = c("LearningSets", "matrix"), 
          definition = function(object, value){
            ##checks??  
            object@learnmatrix <- value
            object
          })

setMethod("method<-", signature = c("LearningSets", "character"), 
          definition = function(object, value){
            ##checks??  
            if(value %in%c("LOOCV", "CV", "MCCV",
                           "bootstrap","custom") ==FALSE )
              stop('Invalid method.')
            object@method <- value
            object
          })

setMethod("ntrain<-", signature = c("LearningSets", "numeric"), 
          definition = function(object, value){
            ##checks??  
            if(value<=0 | length(value!=1 ))
              stop('Invalid "ntrain".')
            object@ntrain <- value
            object
          })

setMethod("iter<-", signature = c("LearningSets", "numeric"), 
          definition = function(object, value){
            ##checks??  
            if(value<=0 | length(value!=1 ))
              stop('Invalid "iter".')
            object@iter <- value
            object
          })

###Constructor
LearningSets<-function(learnmatrix,method,ntrain,iter){
  ###check slots
  if(length(method)!=1 | length(ntrain)!=1 | length(iter)!=1)
    stop("'Method','ntrain' and 'iter' must be of length 1.")
  
  if(nrow(learnmatrix)!=iter)
    stop("Number of rows in 'learnmatrix' must be equal to 'iter'.")
  
  if(any(learnmatrix<0))
    stop("'Learnmatrix' must not contain negative values.")
  
  
  obj<-new("LearningSets",learnmatrix=learnmatrix,method=method,ntrain=ntrain,
           iter=iter)
  return(obj)
}

###############################################################################
###############################################################################
###gene selection class
setClass(Class = "GeneSel", representation(rankings = "list", 
                                           importance = "list", method = "character", 
                                           criterion = "character"))


##accessor(s)
#generic(s)
setGeneric("rankings", function(object) standardGeneric("rankings"))
setGeneric("importance", function(object) standardGeneric("importance"))
setGeneric("method", function(object) standardGeneric("method"))
setGeneric("criterion", function(object) standardGeneric("criterion"))

#method(s)
setMethod("rankings", signature = "GeneSel", definition = 
            function(object) object@rankings)

setMethod("importance", signature = "GeneSel", definition = 
            function(object) object@importance)
setMethod("method", signature = "GeneSel", definition = 
            function(object) object@method)
setMethod("criterion", signature = "GeneSel", definition = 
            function(object) object@criterion)


##replacement function(s)
#generic(s)
setGeneric("rankings<-", function(object, value) standardGeneric("rankings<-"))
setGeneric("importance<-", function(object, value) 
  standardGeneric("importance<-"))
setGeneric("method<-", function(object, value) standardGeneric("method<-"))
setGeneric("criterion<-", function(object, value) 
  standardGeneric("criterion<-"))


#method(s)
setMethod("rankings<-", signature = c("GeneSel", "list"), 
          definition = function(object, value){
            ##checks??  
            object@rankings <- value
            object
          })

setMethod("importance<-", signature = c("GeneSel", "list"), 
          definition = function(object, value){
            ##checks??  
            object@importance <- value
            object
          })

setMethod("method<-", signature = c("GeneSel", "character"), 
          definition = function(object, value){
            ##checks  
            if(value %in% c('fastCox','rowCoxTests','custom'))
              stop('Invalid "method" name.')
            object@method <- value
            object
          })

setMethod("criterion<-", signature = c("GeneSel", "character"), 
          definition = function(object, value){
            ##checks  
            if(value %in% c('pvalue','coefficient','custom'))
              stop('Invalid "criterion".')
            object@criterion <- value
            object
          })


###Constructor
GeneSel<-function(rankings,importance,method,criterion){
  ###checking slots
  if(length(method)!=1 | length(criterion)!=1)
    stop("'Method' and 'criterion' must be of length 1.")
  
  if(length(rankings) != length(importance))
    stop("'Rankings' and 'importance' must be of same length.")
  
  if(any(dim(rankings[[1]])!=dim(importance[[1]])))
    stop("Elements of 'Rankings' and 'importance' must be of same 
         dimension.")
  
  obj<-new("GeneSel",rankings=rankings,importance=importance,method=method,
           criterion=criterion)
  return(obj)
}

###############################################################################
###############################################################################
###filter output class
setClass(Class = "VarSelOut", representation(varsel = "numeric", 
                                             criterion = "character"))

##accessor(s)
#generic(s)
setGeneric("varsel", function(object) standardGeneric("varsel"))
setGeneric("criterion", function(object) standardGeneric("criterion"))

#method(s)
setMethod("varsel", signature = "VarSelOut", definition = 
            function(object) object@varsel)

setMethod("criterion", signature = "VarSelOut", definition = 
            function(object) object@criterion)


##replacement function(s)
#generic(s)
setGeneric("varsel<-", function(object, value) standardGeneric("varsel<-"))
setGeneric("criterion<-", function(object, value) 
  standardGeneric("criterion<-"))

#method(s)
setMethod("varsel<-", signature = c("VarSelOut", "numeric"), 
          definition = function(object, value){
            ##checks??  
            object@varsel <- value
            object
          })


setMethod("criterion<-", signature = c("VarSelOut", "character"), 
          definition = function(object, value){
            ##checks  
            if(value %in% c('pvalue','coefficient','custom'))
              stop('Invalid "criterion".')
            object@criterion <- value
            object
          })




###############################################################################
###############################################################################
###output class of learnsurvival
setClass(Class = "LearnOut", representation(ModelLearnedlist = "list", 
                                            X = "data.frame", y = "Surv", GeneSel = "GeneSel", 
                                            nbgene = "numeric", LearningSets = "LearningSets", 
                                            threshold = "numeric"))


##accessor(s)
#generic(s)
setGeneric("X", function(object) standardGeneric("X"))
setGeneric("y", function(object) standardGeneric("y"))
setGeneric("ModelLearnedlist", function(object) 
  standardGeneric("ModelLearnedlist"))
setGeneric("GeneSel", function(object) standardGeneric("GeneSel"))
setGeneric("nbgene", function(object) standardGeneric("nbgene"))
setGeneric("LearningSets", function(object) standardGeneric("LearningSets"))
setGeneric("threshold", function(object) standardGeneric("threshold"))


#method(s)
setMethod("X", signature = "LearnOut", definition = 
            function(object) object@X)
setMethod("y", signature = "LearnOut", definition = 
            function(object) object@y)
setMethod("ModelLearnedlist", signature = "LearnOut", definition = 
            function(object) object@ModelLearnedlist)
setMethod("GeneSel", signature = "LearnOut", definition = 
            function(object) object@GeneSel)
setMethod("nbgene", signature = "LearnOut", definition = 
            function(object) object@nbgene)
setMethod("LearningSets", signature = "LearnOut", definition = 
            function(object) object@LearningSets)
setMethod("threshold", signature = "LearnOut", definition = 
            function(object) object@threshold)


##replacement function(s)
#generic(s)
setGeneric("X<-", function(object, value) standardGeneric("X<-"))
setGeneric("y<-", function(object, value) standardGeneric("y<-"))
setGeneric("ModelLearnedlist<-", function(object, value) 
  standardGeneric("ModelLearnedlist<-"))
setGeneric("GeneSel<-", function(object, value) 
  standardGeneric("GeneSel<-"))
setGeneric("nbgene<-", function(object, value) standardGeneric("nbgene<-"))
setGeneric("LearningSets<-", function(object, value) 
  standardGeneric("LearningSets<-"))
setGeneric("threshold<-", function(object, value) 
  standardGeneric("threshold<-"))


#method(s)
setMethod("X<-", signature = c("LearnOut", "data.frame"), 
          definition = function(object, value){
            ##checks??  
            object@X <- value
            object
          })

setMethod("y<-", signature = c("LearnOut", "Surv"), 
          definition = function(object, value){
            ##checks??  
            object@y <- value
            object
          })

setMethod("ModelLearnedlist<-", signature = c("LearnOut", "list"), 
          definition = function(object, value){
            ##check
            if(class(value[[1]])!='ModelBase')
              stop('List elements of "ModelLearnedlist" must be of class 
                   "ModelBase" (or a corresponding subclass )')
            object@ModelLearnedlist <- value
            object
          })

setMethod("GeneSel<-", signature = c("LearnOut", "list"), 
          definition = function(object, value){
            ##checks??
            object@GeneSel <- value
            object
          })

setMethod("nbgene<-", signature = c("LearnOut", "numeric"), 
          definition = function(object, value){
            ##checks
            if(value <=0)
              stop('Invalid "nbgene".')
            object@nbgene <- value
            object
          })

setMethod("LearningSets<-", signature = c("LearnOut", "numeric"), 
          definition = function(object, value){
            ##checks->s. class LearningSets
            object@LearningSets <- value
            object
          })

setMethod("threshold<-", signature = c("LearnOut", "numeric"), 
          definition = function(object, value){
            ##checks
            if(value<=0)
              stop('Invalid "threshold".')
            object@threshold <- value
            object
          })


###############################################################################
###############################################################################
### evaluate output class
setClass(Class = "EvaluateOutput", representation(result = "list", 
                                                  measure = "character", additional = "list"))

##accessor(s)
#generic(s)
setGeneric("result", function(object) standardGeneric("result"))
setGeneric("measure", function(object) standardGeneric("measure"))
setGeneric("additional", function(object) standardGeneric("additional"))

#method(s)
setMethod("result", signature = "EvaluateOutput", definition = 
            function(object) object@result)
setMethod("measure", signature = "EvaluateOutput", definition = 
            function(object) object@measure)
setMethod("additional", signature = "EvaluateOutput", definition = 
            function(object) object@additional)


##replacement function(s)
#generic(s)
setGeneric("result<-", function(object, value) standardGeneric("result<-"))
setGeneric("measure<-", function(object, value) standardGeneric("measure<-"))
setGeneric("additional<-", function(object, value)
  standardGeneric("additional<-"))

#method(s)
setMethod("result<-", signature = c("EvaluateOutput", "list"), 
          definition = function(object, value){
            
            ###checks??
            
            object@result <- value
            object
          })

setMethod("measure<-", signature = c("EvaluateOutput", "character"), 
          definition = function(object, value){
            
            ###checks??
            if(value %in% c('UnoC','CvPLogL','HarrellC','PErrC',
                            'custom')==FALSE)
              stop('Invalid measure!')
            object@measure <- value
            object
          })

setMethod("additional<-", signature = c("EvaluateOutput", "list"), 
          definition = function(object, value){
            
            ###checks??
            object@additional <- value
            object
          })



###############################################################################
###############################################################################
### tune output class
setClass(Class = "TuneOut", representation(hypergrid = "data.frame", 
                                           tuneres = "list", method = "character", fold = "numeric", 
                                           X = "data.frame", LearningSets = "LearningSets"))


##accessor(s)
#generic(s)
setGeneric("hypergrid", function(object) standardGeneric("hypergrid"))
setGeneric("tuneres", function(object) standardGeneric("tuneres"))
setGeneric("method", function(object) standardGeneric("method"))
setGeneric("fold", function(object) standardGeneric("fold"))
setGeneric("X", function(object) standardGeneric("X"))
setGeneric("LearningSets", function(object) standardGeneric("LearningSets"))

#method(s)
setMethod("hypergrid", signature = "TuneOut", definition = 
            function(object) object@hypergrid)
setMethod("tuneres", signature = "TuneOut", definition = 
            function(object) object@tuneres)
setMethod("method", signature = "TuneOut", definition = 
            function(object) object@method)
setMethod("fold", signature = "TuneOut", definition = 
            function(object) object@fold)
setMethod("X", signature = "TuneOut", definition = 
            function(object) object@X)
setMethod("LearningSets", signature = "TuneOut", definition = 
            function(object) object@LearningSets)

##replacement function(s)
#generic(s)
setGeneric("hypergrid<-", function(object, value) 
  standardGeneric("hypergrid<-"))
setGeneric("tuneres<-", function(object, value) standardGeneric("tuneres<-"))
setGeneric("method<-", function(object, value) standardGeneric("method<-"))
setGeneric("fold<-", function(object, value) standardGeneric("fold<-"))
setGeneric("X<-", function(object, value) standardGeneric("X<-"))
setGeneric("LearningSets<-", function(object, value) 
  standardGeneric("LearningSets<-"))

#method(s)
setMethod("hypergrid<-", signature = c("TuneOut", "data.frame"), 
          definition = function(object, value){
            ##checks??  
            object@hypergrid <- value
            object
          })

setMethod("tuneres<-", signature = c("TuneOut", "list"), 
          definition = function(object, value){
            ##checks??
            object@tuneres <- value
            object
          })

setMethod("method<-", signature = c("TuneOut", "character"), 
          definition = function(object, value){
            ##checks
            if(value %in%c('plusMinusSurv','coxBoostSurv','customSurv')==FALSE)
              stop('Invalid "method" name.')
            object@method <- value
            object
          })

setMethod("fold<-", signature = c("TuneOut", "numeric"), 
          definition = function(object, value){
            ##checks
            if(value<=0 | length(value)!=1)
              stop('Invalid "fold".')
            object@fold <- value
            object
          })

setMethod("X<-", signature = c("TuneOut", "data.frame"), 
          definition = function(object, value){
            ##checks??
            object@X <- value
            object
          })

setMethod("LearningSets<-", signature = c("TuneOut", "LearningSets"), 
          definition = function(object, value){
            ##checks-> s. class learningsets
            object@LearningSets <- value
            object
          })


###############################################################################
###############################################################################
####evaluation measure classes
##superclass
setClass("MeasureC", representation(type = "character"))


##accessor(s)
#generic(s)
setGeneric("type", function(object) standardGeneric("type"))

#method(s)
setMethod("type", signature = "MeasureC", definition = 
            function(object) object@type)


##replacement function(s)
#generic(s)
setGeneric("type<-", function(object, value) standardGeneric("type<-"))

#method(s)
setMethod("type<-", signature = c("MeasureC", "character"), 
          definition = function(object, value){
            
            if(value %in% c('UnoC','CvPLogL','HarrellC','PErrC',
                            'custom')==FALSE)
              stop('Invalid type!')
            
            object@type <- value
            object
          })

##measurec subclasses (no additional slots!)
setClass("UnoC", representation(measuretype = "character"), 
         prototype(measuretype = "lp"), contains = "MeasureC")
setClass("CvPLogL", representation(measuretype = "character"), 
         prototype(measuretype = "lp"), contains = "MeasureC")

setClass("HarrellC", representation(measuretype = "character"), 
         prototype(measuretype = "lp"), contains = "MeasureC")

setClass("PErrC", representation(measuretype = "character"), 
         prototype(measuretype = "SurvivalProbs"), contains = "MeasureC")



###############################################################################
###############################################################################
###best parameters class 
setClass("BestParameters", representation(param = "ANY", bestind = "numeric", 
                                          perf = "numeric", measure = "MeasureC"))

##accessor(s)
#generic(s)
setGeneric("param", function(object) standardGeneric("param"))
setGeneric("bestind", function(object) standardGeneric("bestind"))
setGeneric("perf", function(object) standardGeneric("perf"))
setGeneric("measure", function(object) standardGeneric("measure"))

#method(s)
setMethod("param", signature = "BestParameters", definition = 
            function(object) object@param)
setMethod("bestind", signature = "BestParameters", definition = 
            function(object) object@bestind)
setMethod("perf", signature = "BestParameters", definition = 
            function(object) object@perf)
setMethod("measure", signature = "BestParameters", definition = 
            function(object) object@measure)


##replacement function(s)
#generic(s)
setGeneric("param<-", function(object, value) standardGeneric("param<-"))
setGeneric("bestind<-", function(object, value) standardGeneric("bestind<-"))
setGeneric("perf<-", function(object, value) standardGeneric("perf<-"))
setGeneric("measure<-", function(object, value) standardGeneric("measure<-"))


#method(s)
setMethod("param<-", signature = c("BestParameters", "ANY"), 
          definition = function(object, value){
            ##checks??  
            
            object@param <- value
            object
          })

setMethod("bestind<-", signature = c("BestParameters", "numeric"), 
          definition = function(object, value){
            ##checks  
            if(length(value)!=1)
              stop('"bestind" must be of length 1.')
            
            object@bestind <- value
            object
          })

setMethod("perf<-", signature = c("BestParameters", "numeric"), 
          definition = function(object, value){
            ##checks??  
            object@perf <- value
            object
          })

setMethod("measure<-", signature = c("BestParameters", "MeasureC"), 
          definition = function(object, value){
            ##checks->  s. MeasureC
            object@measure <- value
            object
          })


###############################################################################
###############################################################################
###compare measures
##superclass
setClass(Class="MeasureOut",representation(estimate="numeric",
                                           conf.int="numeric"))


##accessor(s)
#generic(s)
setGeneric("estimate", function(object) standardGeneric("estimate"))
setGeneric("conf.int", function(object) standardGeneric("conf.int"))

#method(s)
setMethod("estimate", signature = "MeasureOut", definition = 
            function(object) object@estimate)
setMethod("conf.int", signature = "MeasureOut", definition = 
            function(object) object@conf.int)


##replacement function(s)
#generic(s)
setGeneric("estimate<-", function(object, value) standardGeneric("estimate<-"))
setGeneric("conf.int<-", function(object, value) standardGeneric("conf.int<-"))

#method(s)
setMethod("estimate<-", signature = c("MeasureOut", "numeric"), 
          definition = function(object, value){
            ##checks??  
            if(length(value)!=1)
              stop('"estimate" must be of length 1.')
            
            object@estimate <- value
            object
          })

setMethod("conf.int<-", signature = c("MeasureOut", "numeric"), 
          definition = function(object, value){
            ##checks??  
            if(length(value)!=2)
              stop('"conf.int" must be of length 2.')
            
            object@conf.int <- value
            object
          })






##sub-classes (no additional slots)
###class NRI
setClass(Class="IDI",contains="MeasureOut")

####class ID
setClass(Class="NRI",contains="MeasureOut")

####class MIRS
setClass(Class="MIRS",contains="MeasureOut")

#####class CDelta
setClass(Class="CDelta",contains="MeasureOut")