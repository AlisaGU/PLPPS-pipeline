# Title: Various gene selection methods.
#
# Author: M. Slawski, adapted from A-L Boulesteix
# Email: <Martin.Slawski@campus.lmu.de>
# Date of creation: 26.9.2007
#
# Brief description:
#   Returns an object of class 'GeneSel'.
#
# Arguments:
#   -X: matrix of variables (rows are observations,columns are variables)
#   -y: survival response of class Surv 
#   -f: formula for survival model
#   -LearningSets: output-object of method GenerateLearningsets
#	(splits of dataset into training/test sets)
#   -method: variable selection method to be used
#   -trace: print some additional information 
#   -criterion: pvalue or coefficients, which shall be returned by 
#	the filter method?
#
# Value:
#   GeneSel
#
###############################################################################

setGeneric("geneSelection", function(X, y, ...) 
  standardGeneric("geneSelection"))

setMethod("geneSelection", signature(X = "data.frame", y = "Surv"), 
          function(X, y, LearningSets, method = c("fastCox"), 
                   criterion = c("coefficient"), trace = TRUE, ...) {
            nrx <- nrow(X)
            if (nrx != nrow(y)) 
              stop("Number of rows of 'X' must agree with length of y \n")
            tempcall <- as.character(match.call())
            
            if (missing(LearningSets)) {
              warning("Argument 'LearningSets' is missing; set to a row vector
                      with entries '1:nrow(X)' \n")
              learnmatrix <- matrix(1:nrx, ncol = nrx)
            } else {
              learnmatrix <- LearningSets@learnmatrix
              if (ncol(learnmatrix) > nrx) 
                stop("'LearningSets' do not match the input data \n")
            }
            
            niter <- nrow(learnmatrix)
            p <- ncol(X)
            outrankings <- outimportance <- matrix(nrow = niter, ncol = p)
            rankings <- importance <- matrix(0, niter, p)
            selfun <- switch(method, fastCox = fastCox, 
                             stop("Invalid 'method' specified\n"))
            
            for (i in 1:niter) {
              if (trace) 
                cat("geneSelection: iteration", i, "\n")
              
              outporg <- selfun(X, y, learnind = learnmatrix[i, ], 
                                criterion = criterion, ...)
              outp <- outporg@varsel
              decr <- outporg@criterion != "pvalue"
              outrankings[i, ] <- ord <- order(outp, decreasing = decr)
              outimportance[i, ] <- outp[ord]
            }
            
            colnames(outrankings) <- paste("rank", 1:p, sep = "")
            colnames(outimportance) <- paste("gene", ord, sep = "")
            rownames(outrankings) <- rownames(outimportance) <- 
              paste("iter.", 1:niter, sep = "")
            rankings <- importance <- list()
            rankings[[1]] <- outrankings
            importance[[1]] <- outimportance
            
            new("GeneSel", rankings = rankings, importance = importance, 
                method = method, criterion = criterion)
          })

setMethod("geneSelection", signature(X = "ExpressionSet", y = "character"), 
          function(X, y, ... ) {
            Xdat <- as.data.frame(t(exprs(X)), check.names = FALSE)
            geneSelection(X = Xdat, y = .fetchyFromEset(X,y), ...)
          })

setMethod("geneSelection", signature(X = "ExpressionSet", y = "Surv"), 
          function(X, y,...) {
            geneSelection(X = as.data.frame(t(exprs(X)), check.names = FALSE), 
                          y = y, ... )}) 

setMethod("geneSelection", signature(X = "matrix", y = "Surv"), function(X, y,
                                                                         ...) {
  geneSelection(X = as.data.frame(X, check.names = FALSE), 
                y = y, ... )
}) 