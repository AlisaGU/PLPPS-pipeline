# Author: M. Slawski, modified by Levi Waldron
# Email: <Martin.Slawski@campus.lmu.de>
# Date of creation: 9.10.2007
# Adapted 01.16.2012 for survival
#
# Brief description:
#   Returns a list of objects of class "ModelLearned", for various 
#   LearningSets. "General" function of the package.
#
# Further comments and notes:
#   evaluationSurv.r, geneSelection.r learnSurvival.r
#
# Arguments
#   -X: matrix of variables (rows are observations,columns are variables)
#   -y: survival response of class Surv 
#   -f: formula for survival model
#   -LearningSets: output-object of method GenerateLearningsets (splits of 
#	dataset into training/test sets)
#   -GeneSel: return object of method "GeneSelection"
#   -nbgene: number of genes to be selected
#   -survmethod: survival method to be used (simple "character")
#   -tuneres: return-object of method "tune"
#   -trace: print some additional information
#   -innercv: logical indicating whether learnSurvival is used in an inner 
#	loop (X will not be stored)
#   -addtune: list with additional parameters for evaluating tuneres
#   -threshold: threshold for pvalue or coefficient when selecting genes
#
# Value:
#   LearnOut
#
###############################################################################

setGeneric("learnSurvival", function(X, y, ...) 
  standardGeneric("learnSurvival"))

setMethod("learnSurvival", signature(X = "data.frame", y = "Surv"), 
          function(X, y, LearningSets, GeneSel, GeneSellist = list(), nbgene, 
                   survmethod, tuneres, tuninglist = list(), trace = TRUE, 
                   innercv = FALSE, addtune = list(), threshold = NULL, ...) {
            if (!missing(nbgene) & !is.null(threshold)) 
              stop("You can either specify a fixed number of selected genes
                   (nbgene) or a threshold for selecting genes.")
            
            ll <- list(...)
            if (missing(survmethod)) 
              stop("argument 'survmethod' is missing \n")
            if (missing(LearningSets)) {
              # warning("Argument 'LearningSets' is missing; set to a row
              #vector with entries '1:nrow(X)' \n")
              learnmatrix <- matrix(1:nrow(X), nrow = 1)
              LearningSets <- new("LearningSets", learnmatrix =
                                    learnmatrix, method = "usealldata", 
                                  ntrain = ncol(learnmatrix), iter = 1)
            } else {
              learnmatrix <- LearningSets@learnmatrix
              if (ncol(learnmatrix) > nrow(X)) 
                stop("'LearningSets' do not match the input data \n")
            }
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
                stop("'nbgene' greater than the number of all genes \n")
            } else nbgene <- ncol(X)
            if (missing(tuneres)) {
              if (!missing(tuninglist) && length(tuninglist) != 0) {
                tuninglist$X <- X
                tuninglist$y <- y
                tuninglist$survmethod <- survmethod
                if (!missing(LearningSets)) 
                  tuninglist$LearningSets <- LearningSets
                if (!missing(GeneSel)) {
                  tuninglist$GeneSel <- GeneSel
                  tuninglist$nbgene <- nbgene
                }
                if (!is.list(tuninglist$grid)) 
                  stop("Invalid specification of 'tuninglist'. Grid must
                       itself be a list \n")
                tuneres <- do.call("tune", args = c(tuninglist, ll))
              }
            }
            if (!missing(tuneres)) {
              if (length(tuneres@tuneres) != nrow(learnmatrix)) 
                stop("object 'tuneres' does not match 'LearningSets' \n")
              if (!grep(tuneres@method, survmethod, ignore.case = TRUE)) 
                stop("object 'tuneres' does not match the 
                     chosen survmethod. \n")
              besthyperpar <- list()
              addtune <- c(addtune, list(...))
              addtune$tuneres <- tuneres
              for (indloop in 1:length(tuneres@tuneres)) {
                addtune$res.ind <- indloop
                besthyperpar[[indloop]] <- 
                  do.call("getBestParameters", args = addtune)@param
              }
              
            }
            survoutlist <- vector(mode = "list", length = nrow(learnmatrix))
            if (missing(GeneSel)) {
              if (missing(tuneres)) {
                for (i in 1:nrow(learnmatrix)) {
                  if (trace) 
                    cat("iteration", i, "\n")
                  survoutlist[[i]] <- do.call(survmethod, args = 
                                                c(list(X = X, y = y, 
                                                       learnind = learnmatrix[i,]),ll))
                }
              } else {
                arglist <- expression(c(ll, list(X = X, y = y, 
                                                 learnind = learnmatrix[i, 
                                                                        ]), besthyperpar[[i]]))
                funtocall <- survmethod
                for (i in 1:nrow(learnmatrix)) {
                  if (trace) 
                    cat("iteration", i, "\n")
                  survoutlist[[i]] <- do.call(funtocall, 
                                              args = eval(arglist))
                }
              }
            } else {
              ranks <- GeneSel@rankings
              imps <- GeneSel@importance
              criterion <- GeneSel@criterion
              if (missing(tuneres)) {
                arglist <- expression(c(ll, list(X = Xi, y = y, 
                                                 learnind = learnmatrix[i, 
                                                                        ])))
                funtocall <- survmethod
              } else {
                arglist <- expression(c(ll, list(X = Xi, y = y, 
                                                 learnind = learnmatrix[i, 
                                                                        ]), besthyperpar[[i]]))
                funtocall <- survmethod
              }
              if (!is.null(threshold)) {
                if (length(ranks) > 1) {
                  for (i in 1:nrow(learnmatrix)) {
                    if (trace) 
                      cat("iteration", i, "\n")
                    seli <- c()
                    for (j in 1:length(ranks)) {
                      rankj <- ranks[[j]][i, ]
                      impj <- imps[[j]][i, ]
                      if (criterion == "coefficient") 
                        impj <- impj[impj > threshold] 
                      else impj <- impj[impj < threshold]
                      
                      nbgene2 <- min(length(impj), nbgene)
                      seli <- c(seli, rankj[1:nbgene2])
                    }
                    Xi <- X[, seli, drop = FALSE]
                    survoutlist[[i]] <- do.call(funtocall, 
                                                args = eval(arglist))
                  }
                } else {
                  ranks <- ranks[[1]]
                  imps <- imps[[1]]
                  for (i in 1:nrow(learnmatrix)) {
                    if (trace) 
                      cat("iteration", i, "\n")
                    impi <- imps[i, ]
                    
                    if (criterion == "coefficient") 
                      impi <- impi[impi > threshold] 
                    else impi <- impi[impi < threshold]
                    
                    
                    nbgene2 <- min(length(impi), nbgene)
                    seli <- ranks[i, 1:nbgene2]
                    Xi <- X[, seli, drop = FALSE]
                    survoutlist[[i]] <- 
                      do.call(funtocall, args = eval(arglist))
                  }
                }
              } else {
                if (length(ranks) > 1) {
                  for (i in 1:nrow(learnmatrix)) {
                    if (trace) 
                      cat("iteration", i, "\n")
                    seli <- c()
                    for (j in 1:length(ranks)) {
                      rankj <- ranks[[j]][i, ]
                      seli <- c(seli, rankj[1:nbgene])
                    }
                    seli <- unique(seli)
                    Xi <- X[, seli, drop = FALSE]
                    survoutlist[[i]] <- 
                      do.call(funtocall, args = eval(arglist))
                  }
                } else {
                  ranks <- ranks[[1]]
                  imps <- imps[[1]]
                  for (i in 1:nrow(learnmatrix)) {
                    if (trace) 
                      cat("iteration", i, "\n")
                    seli <- ranks[i, 1:nbgene]
                    Xi <- X[, seli, drop = FALSE]
                    survoutlist[[i]] <- 
                      do.call(funtocall, args = eval(arglist))
                  }
                }
              }
            }
            
            if (class(try(checkGeneSel <- is.null(GeneSel), silent = TRUE)) == 
                "try-error") 
              GeneSel <- new("GeneSel")
            if (is.null(threshold)) 
              threshold = -1
            ### return a LearnOut object
            if (!innercv) 
              new("LearnOut", ModelLearnedlist = survoutlist, X = X, 
                  y = y, nbgene = nbgene, GeneSel = GeneSel, 
                  LearningSets = LearningSets, threshold = threshold) 
            else new("LearnOut", ModelLearnedlist = survoutlist, 
                     X = as.data.frame(0), y = y, nbgene = nbgene, 
                     GeneSel = GeneSel, LearningSets = LearningSets, 
                     threshold = threshold)
            })

setMethod("learnSurvival", signature(X = "ExpressionSet", y = "character"),
          function(X, y, ...) {
            Xdat <- as.data.frame(t(exprs(X)), check.names = FALSE)
            learnSurvival(X = Xdat, y = .fetchyFromEset(X,y), ...)
          })

setMethod("learnSurvival", signature(X = "ExpressionSet", y = "Surv"),
          function(X, y, ...) { 
            learnSurvival(X = as.data.frame(t(exprs(X)), check.names = FALSE), 
                          y = y, ...)}) 

setMethod("learnSurvival", signature(X = "matrix", y = "Surv"), 
          function(X, y, ...) { 
            learnSurvival(X = as.data.frame(X, check.names = FALSE), 
                          y = y, ...)
          }) 
