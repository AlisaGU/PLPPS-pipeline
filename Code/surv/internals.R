.fetchyFromEset <- function(X, y) {
  if (!(length(y) %in% 1:2) || any(y %in% varLabels(X)) == FALSE) 
    stop("if X is ExpressionSet and y is character, y must have length 
         one or two and all its elements must correspond to one 
         of varLabels(X).")
  if (length(y) == 1) 
    ysurv <- X[[y]] else {
      if (class(X[[y[2]]]) == "factor") 
        stop("Survival indicator is of class \"factor\" and the
             exact coding is unclear. Please provide a logical 
             or integer variable.")
      ysurv <- Surv(X[[y[1]]], X[[y[2]]])
    }
  ysurv
}