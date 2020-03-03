setMethod("show", signature = "LearningSets", function(object) {
  cat("learningset mode: ", object@method, "\n")
  cat("number of LearningSets: ", object@iter, "\n")
  cat("(maximum) number of observations per learning set: ",
      object@ntrain, "\n")
})

setMethod("show", signature = "GeneSel", function(object) {
  cat("gene selection performed with '", 
      object@method, "'", "\n", sep = "")
  cat("criterion used :'", object@criterion, "'", "\n", sep = "")
  if (length(object@rankings) == 1) 
    ngenes <- ncol(object@rankings[[1]]) 
  else ngenes <- unlist(lapply(object@rankings, ncol))[1]
  cat("number of genes: ", ngenes, "\n")
  if (length(object@rankings) == 1) 
    niter <- nrow(object@rankings[[1]]) 
  else niter <- unlist(lapply(object@rankings, nrow))[1]
  cat("number of different LearningSets: ", niter, "\n")
})

setMethod("show", signature = "TuneOut", function(object) {
  cat("tuning for '", object@method, "'", "\n", sep = "")
  cat("hyperparameters tuned: \n")
  cat(colnames(object@hypergrid), sep = ",")
  cat("\n")
  cat("CV folds used:", object@fold, "\n")
  cat("\n")
})
