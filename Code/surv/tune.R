### filename: tune.r
### Title: Function to tune different classication methods.
###
###
### Author: C.Bernau, adapted from M. Slawski
### email: <bernau@ibe.med.uni-muenchen.de>
### date of creation: 24.4.2012
#
### Brief description:
#
#   Returns an object of class "TuneOut",
#   based on results of "evaluationsurv".
#
### Further comments and notes:
#   s. evaluationSurv.r
#   s. geneSelection.r
#   s. learnSurvival.r
#
##########################
#arguments
#-X: matrix of variables (rows are observations,columns are variables)
#-y: survival response of class Surv 
#-f: formula for survival model
#-LearningSets: output-object of method generateLearningsets
# (splits of dataset into training/test sets)
#-GeneSel: return object of method "geneSelection"
#-nbgene: number of genes to be selected
#-survmethod: survival method to be used (simple "character")
#-fold: number of cross-validation folds
#-strat: use stratified cv
#-grid: list of grids for the tuning parameter
#-trace: print some additional information
#-threshold: threshold for pvalue or coefficient when selecting genes
#-strat: use stratified cv
#
####################
#return-class
#TuneOutws
#
###**************************************************************************###
### generic

setGeneric("tune", function(X, y, ...) standardGeneric("tune"))

### X=data.frame, y=Surv, LearningSets=LearningSets, measure=missing (no direct
### evaluation)

setMethod("tune", signature(X = "data.frame", y = "Surv"), 
          function(X, y, LearningSets, measure=NULL, GeneSel, 
                   GeneSellist = list(), nbgene, survmethod, fold = 3, 
                   strat = FALSE, grids = list(), trace = TRUE, 
                   threshold = NULL, ...) {
            
            
            if (is.null(LearningSets)) LearningSets <- new("LearningSets", 
                                                           learnmatrix = matrix(1:nrow(X), nrow = 1), 
                                                           method = "whole data set", 
                                                           ntrain = nrow(X), iter = 1)
            
            if(class(LearningSets)!="LearningSets") stop("Parameter 
                                                         'LearningSets' must be of class 'LearningSets'")
            
            if(is.null(measure)) directeval<-FALSE
            
            else directeval<-TRUE
            
            if (!missing(nbgene) & is.null(threshold) == FALSE) 
              stop("You can either specify a fixed number of 
                   selected genes (nbgene) or a threshold for 
                   selecting genes.")
            
            addtune <- list(...)
            
            if(is.null(addtune$tuneres)==FALSE) {
              warning('Argument "tuneres" will be ignored.')
              addtune$tuneres<-NULL}
            
            learnmatrix <- LearningSets@learnmatrix
            if (ncol(learnmatrix) > nrow(X)) 
              stop("'LearningSets' do not match the input data \n")
            
            if (missing(GeneSel)) {
              if (!missing(GeneSellist) && length(GeneSellist) != 0) {
                GeneSellist$X <- X
                GeneSellist$y <- y
                if (!missing(LearningSets)) 
                  GeneSellist$LearningSets <- LearningSets
                GeneSel <- do.call("geneSelection", args = GeneSellist)
              }
            } else {
              if (class(GeneSel) != "GeneSel") 
                stop("'GeneSel' must be of class 'GeneSel' \n")
              ngenes <- ncol(GeneSel@rankings[[1]])
              niterGeneSel <- nrow(GeneSel@rankings[[1]])
              if (ngenes != ncol(X)) 
                stop("object 'GeneSel' does not match the input data \n")
              if (niterGeneSel != nrow(learnmatrix)) 
                stop("object 'GeneSel' does not match 'LearningSets' \n")
            }
            
            if (!missing(nbgene)) {
              if (nbgene > ncol(X)) 
                stop("'nbgene' greater than the number all genes \n")
              if (GeneSel@method == "lasso" | GeneSel@method == "elasticnet")
              {
                if (min(unlist(lapply(GeneSel@importance, 
                                      function(y) apply(y, 1, function(x) 
                                        which(x == 0)))))) {
                  warning("'nbgene' greater than number of 
                          nonzero-coefficients in lasso/elasticnet
                          on at least one training set.")
                }
                }
                } else nbgene <- ncol(X)
            
            ll <- list(...)
            
            if(is.null(ll$tuneres)==FALSE) {
              warning('Argument "tuneres" will be ignored.')
              ll$tuneres<-NULL}
            
            ll$innercv <- TRUE
            survname <- survmethod
            print(survname)
            if (length(grids) == 0) {
              grids <- switch(survname, penalizedSurv = list(lambda = 
                                                               c(0.01, 0.1, 1, 10, 100, 1000, 10000)),
                              coxBoostSurv = list(stepno = c(10 * (1:10))))
            }
            
            if (length(grids) == 0) 
              stop("'survival method' does not need any tuning \n")
            innerlength <- unlist(lapply(grids, length))
            if (any(innerlength == 0)) 
              stop("Invalid grids specified \n")
            
            ## sorting grids
            grids <- lapply(grids, sort)
            hypergrid <- expand.grid(grids)
            
            tunereslist <- vector(mode = "list", length = nrow(learnmatrix))
            
            if (missing(GeneSel)) {
              tuneresi <- list()
              for (i in 1:nrow(learnmatrix)) {
                if (trace) 
                  cat("tuning iteration", i, "\n")
                Xi <- X[learnmatrix[i, ], , drop = FALSE]
                yi <- y[learnmatrix[i, ], ]
                lsi <- generateLearningsets(y = yi, method = "CV", 
                                            fold = fold, strat = strat)
                
                for (k in 1:nrow(hypergrid)) {
                  lsurvk <- do.call("learnSurvival", 
                                    args = c(list(X = Xi, y = yi, 
                                                  LearningSets = lsi, trace = FALSE, 
                                                  survmethod = survmethod), 
                                             as.list(data.frame(hypergrid[k, , drop = 
                                                                            FALSE])), ll))
                  
                  if(directeval){
                    evalklist <- addtune
                    evalklist$X <- Xi
                    evalklist$measure <- measure
                    evalklist$object <- lsurvk
                    evalk <- do.call("evaluate", args = evalklist)
                    tuneresi[[k]] <- evalk}
                  
                  else tuneresi[[k]] <- lsurvk
                }
                tunereslist[[i]] <- tuneresi
              }
            } else {
              
              ranks <- GeneSel@rankings
              imps <- GeneSel@importance
              criterion <- GeneSel@criterion
              if (is.null(threshold) == FALSE) {
                if (length(ranks) > 1) {
                  for (i in 1:nrow(learnmatrix)) {
                    tuneresi <- list()
                    if (trace) 
                      cat("tuning iteration", i, "\n")
                    seli <- c()
                    for (j in 1:length(ranks)) {
                      rankj <- ranks[[j]][i, ]
                      impj <- imps[[j]][i, ]
                      
                      if (criterion == "coefficient") 
                        impj <- impj[impj > threshold] else impj <- 
                        impj[impj < threshold]
                      
                      nbgene <- min(length(impj), nbgene)
                      seli <- c(seli, rankj[1:nbgene])
                    }
                    seli <- unique(seli)
                    Xi <- X[learnmatrix[i, ], seli, drop = FALSE]
                    yi <- y[learnmatrix[i, ], ]
                    lsi <- generateLearningsets(y = yi, method = "CV",
                                                fold = fold, strat = strat)
                    for (k in 1:nrow(hypergrid)) {
                      lsurvk <- do.call("learnSurvival", args =
                                          c(list(X = Xi, y = yi, 
                                                 LearningSets = lsi, trace = FALSE, 
                                                 survmethod = survmethod), 
                                            as.list(data.frame(hypergrid[k, , 
                                                                         drop = FALSE])), ll))
                      
                      if(directeval){
                        evalklist <- addtune
                        evalklist$X <- Xi
                        evalklist$measure <- measure
                        evalklist$object <- lsurvk
                        evalk <- do.call("evaluate", args = 
                                           evalklist)
                        tuneresi[[k]] <- evalk}
                      
                      else tuneresi[[k]] <- lsurvk
                    }
                    tunereslist[[i]] <- tuneresi
                  }
                } else {
                  ranks <- ranks[[1]]
                  imps <- imps[[1]]
                  for (i in 1:nrow(learnmatrix)) {
                    tuneresi <- list()
                    if (trace) 
                      cat("tuning iteration", i, "\n")
                    impi <- imps[i, ]
                    
                    if (criterion == "coefficient") 
                      impi <- impi[impi > threshold] else impi <- 
                      impi[impi < threshold]
                    
                    nbgene <- min(length(impi), nbgene)
                    seli <- ranks[i, 1:nbgene]
                    Xi <- X[learnmatrix[i, ], seli, drop = FALSE]
                    yi <- y[learnmatrix[i, ]]
                    lsi <- generateLearningsets(y = yi, method = "CV",
                                                fold = fold, 
                                                strat = strat)
                    for (k in 1:nrow(hypergrid)) {
                      lsurvk <- do.call("learnSurvival", args = 
                                          c(list(X = Xi, y = yi, 
                                                 LearningSets = lsi, trace = FALSE,
                                                 survmethod = survmethod), 
                                            as.list(data.frame(hypergrid[k, , 
                                                                         drop = FALSE])), ll))
                      
                      if(directeval){
                        evalklist <- addtune
                        evalklist$X <- Xi
                        evalklist$measure <- measure
                        evalklist$object <- lsurvk
                        evalk <- do.call("evaluate", 
                                         args = evalklist)
                        tuneresi[[k]] <- evalk}
                      
                      else tuneresi[k] <- lsurvk
                    }
                    
                    
                    tunereslist[[i]] <- tuneresi
                  }
                }
              } else {
                if (length(ranks) > 1) {
                  for (i in 1:nrow(learnmatrix)) {
                    tuneresi <- list()
                    if (trace) 
                      cat("tuning iteration", i, "\n")
                    seli <- c()
                    for (j in 1:length(ranks)) {
                      rankj <- ranks[[j]][i, ]
                      seli <- c(seli, rankj[1:nbgene])
                    }
                    seli <- unique(seli)
                    Xi <- X[learnmatrix[i, ], seli, drop = FALSE]
                    yi <- y[learnmatrix[i, ]]
                    lsi <- generateLearningsets(y = yi, method = "CV", 
                                                fold = fold, 
                                                strat = strat)
                    perf <- double(nrow(hypergrid))
                    for (k in 1:nrow(hypergrid)) {
                      lsurvk <- do.call("learnSurvival", 
                                        args = c(list(X = Xi, y = yi, 
                                                      LearningSets = lsi, trace = FALSE, 
                                                      survmethod = survmethod), 
                                                 as.list(data.frame(hypergrid[k, ,
                                                                              drop = FALSE])), ll))
                      
                      if(directeval){
                        evalklist <- addtune
                        evalklist$X <- Xi
                        evalklist$measure <- measure
                        evalklist$object <- lsurvk
                        evalk <- do.call("evaluate", 
                                         args = evalklist)
                        tuneresi[[k]] <- evalk}
                      
                      
                      else tuneresi[k] <- lsurvk
                    }
                    tunereslist[[i]] <- tuneresi
                  }
                } else {
                  ranks <- ranks[[1]]
                  for (i in 1:nrow(learnmatrix)) {
                    tuneresi <- list()
                    if (trace) 
                      cat("tuning iteration", i, "\n")
                    seli <- ranks[i, 1:nbgene]
                    Xi <- X[learnmatrix[i, ], seli, drop = FALSE]
                    yi <- y[learnmatrix[i, ]]
                    lsi <- generateLearningsets(y = yi, method = "CV",
                                                fold = fold, 
                                                strat = strat)
                    perf <- double(nrow(hypergrid))
                    for (k in 1:nrow(hypergrid)) {
                      lsurvk <- do.call("learnSurvival", 
                                        args = c(list(X = Xi, y = yi, 
                                                      LearningSets = lsi, trace = FALSE, 
                                                      survmethod = survmethod), 
                                                 as.list(data.frame(hypergrid[k, , 
                                                                              drop = FALSE])), ll))
                      
                      if(directeval){
                        evalklist <- addtune
                        evalklist$X <- Xi
                        evalklist$measure <- measure
                        evalklist$object <- lsurvk
                        evalk <- do.call("evaluate", 
                                         args = evalklist)
                        tuneresi[[k]] <- evalk}
                      
                      else tuneresi[k] <- lsurvk
                    }
                    tunereslist[[i]] <- tuneresi
                  }
                }
              }
            }
            
            return(new("TuneOut", hypergrid = hypergrid, tuneres = tunereslist, 
                       method = survmethod, 
                       fold = fold, X = X, LearningSets = LearningSets))
            
              })



### X=ExpressionSet, y=Surv, LearningSets=ANY

setMethod("tune", signature(X = "ExpressionSet", y = "Surv"), function(X, y, 
                                                                       LearningSets = NULL, measure = NULL, ...) {
  Xdat <- as.data.frame(t(exprs(X)), check.names = FALSE)
  tune(X = Xdat, y = y, LearningSets = LearningSets, 
       measure = measure,...)
})

setMethod("tune", signature(X = "ExpressionSet", y = "character"), 
          function(X, y, LearningSets = NULL, measure = NULL, ...) { 
            Xdat <- as.data.frame(t(exprs(X)), check.names = FALSE)
            if (!(length(y) %in% 1:2) || any(y %in% varLabels(X)) == FALSE) 
              stop("if X is ExpressionSet and y is character, y must have 
                   length one or two and all its elements must 
                   correspond to one of varLabels(X).")
            
            if (!(length(y) %in% 1:2) || any(y %in% varLabels(X)) == FALSE) 
              stop("if X is ExpressionSet and y is character, y must have 
                   length one or two and all its elements must 
                   correspond to one of varLabels(X).")
            
            if (length(y) == 1) 
              ysurv <- X[[y]] else {
                if (class(X[[y[2]]]) == "factor") 
                  stop("Survival indicator is of class \"factor\" 
                       and the exact coding is unclear. 
                       Please provide
                       a logical or integer variable.")
                ysurv <- Surv(X[[y[1]]], X[[y[2]]])
              }
            tune(X = Xdat, y = ysurv, LearningSets = LearningSets, 
                 measure = measure, ...) 
          })

### X=matrix, y=Surv, LearningSets=ANY

setMethod("tune", signature(X = "matrix", y = "Surv"), 
          function(X, y, LearningSets = NULL, measure = NULL, ...) {
            tune(X = as.data.frame(X, check.names = FALSE), y = y, 
                 LearningSets = LearningSets, 
                 measure = measure, ...)
          })

