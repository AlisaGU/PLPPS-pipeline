# Filename: GenerateLearningsets.r
# Title: Function to prepare different LearningSets, e.g.
#        - LOOCV
#        - MCCV
#        - generateBootstrap
#
#
# Author: M. Slawski, adapted from A-L Boulesteix
# Email: <Martin.Slawski@campus.lmu.de>
# Date of creation: 19.9.2007
#
# Brief description:
#   Returns an object of class LearningSets.
#
# Further comments and notes:
#   BoostrapCV is, though highly interesting, not included because it
#   does not fit into the general design.
#
# Arguments:
#   -n: number of observations
#   -y: survival response (class Surv)
#   -fold: number of folds in CV
#   -niter: number of iterations
#   -ntrain: number of training observations
#   -strat: use stratified resampling 
#   
# Value:
#   LearningSets
#
###############################################################################

generateLearningsets <- function(n, y, method = c("LOOCV", "CV", "MCCV",
                                                  "bootstrap"), fold = NULL, niter = NULL, ntrain = NULL, 
                                 strat = FALSE) {
  
  ## if(class(y)=='Surv') ## LW: got rid of this because y can be missing, and
  ## because y[,2] (censoring) is used below.  y<-y[,1]
  
  if (!missing(n)) {
    if (length(n) != 1 || n < 0) 
      stop("'n' must be a positive integer ! \n")
    n <- as.integer(n)
    if (!is.null(fold) && n <= fold) 
      stop("'n' is too small \n")
    if (!is.null(ntrain) && n <= ntrain) 
      stop("'n' is too small \n")
  }
  
  ## For class(y) == 'Surv', convert y to a factor of the censoring variable.
  ## The effect of this is to use the censoring variable only for 
  ##stratification.
  if (!missing(y)) {
    if (identical(all.equal(class(y), "Surv"), TRUE)) {
      y <- y[, 2]  #censoring variable
      y <- factor(y)
      #if (strat) 
      #   warning("For survival response, stratification is done by 
      #censoring only.")
    }
  }
  
  if (missing(n) & missing(y)) 
    stop("At least one of 'n' or 'y' mus be given \n")
  if (!missing(y)) 
    n <- length(y)
  
  if (strat & missing(y)) 
    stop("If 'strat=TRUE', 'y' (class memberships) must be given \n")
  
  method <- match.arg(method) 
  
  if (method == "MCCV") {
    if (is.null(niter) | is.null(ntrain)) 
      stop("With the MCCV method, arguments niter and ntrain 
           should be given.")
    
    if(!missing(y)){
      taby <- table(y)
      prop <- taby/sum(taby)
      classize <- .roundVector(prop * ntrain, ntrain)
      if (any(classize < 1) & strat==TRUE) {
        warning("Not enough censored or uncensored observations for 
                stratification. Argument 'strat' set to 
                'FALSE'. \n")
        srat<-FALSE}
    }
    
    if (strat) {
      
      indlist <- sapply(names(taby), function(z) which(y == z), 
                        simplify = FALSE)
      learnmatrix <- matrix(nrow = niter, ncol = ntrain)
      lower <- cumsum(c(1, classize[-length(classize)]))
      upper <- cumsum(classize)
      for (i in 1:length(indlist)) learnmatrix[, lower[i]:upper[i]] <-
        t(replicate(niter, 
                    sample(indlist[[i]], classize[i], 
                           replace = FALSE)))
    } else learnmatrix <- t(replicate(niter, sample(n, ntrain, 
                                                    replace = FALSE)))
    
    } else if (method == "CV") {
      if (is.null(niter)) 
        niter <- 1
      if (is.null(fold)) 
        stop("With the CV method, argument 'fold' must be given.")
      
      if(!missing(y)){
        taby <- table(y)
        prop <- taby/sum(taby)
        siz <- n - floor(n/fold)
        classize <- .roundVector(prop * siz, siz)
        if (any(taby < fold) & strat==TRUE) {
          warning("Not enough censored or uncensored observations for 
                  stratification. Argument 'strat' set to 
                  'FALSE'. \n")
          srat<-FALSE}
      }
      
      
      if (!strat) {
        if (fold == n) 
          method <- "LOOCV" else {
            size <- n/fold
            # if (size < 5) stop('argument 'fold' is too large;
            #The ratio of no.
            # observations/fold should be > 5. \n')
            learnmatrix <- matrix(0, niter * fold, n - floor(size))
            size.int <- floor(size)
            size.vector <- rep(size.int, fold)
            
            if (size.int != size) 
              size.vector[1:((size - size.int) * fold)] <- 
              size.vector[1:((size - 
                                size.int) * fold)] + 1
            
            group.index <- c()
            for (j in 1:fold) group.index <- c(group.index, 
                                               rep(j, size.vector[j]))
            
            for (i in 1:niter) {
              group.index <- group.index[sample(n, n, replace = FALSE)]
              for (j in 1:fold) {
                whichj <- which(group.index == j)
                learnmatrix[j + (i - 1) * fold, 1:length(whichj)] <-
                  whichj
              }
            }
            learnmatrix <- learnmatrix[, 1:max(size.vector), drop = FALSE]
            
            if (size.int != size) 
              learnmatrix <- t(apply(learnmatrix, 1, 
                                     function(z) setdiff(0:n, z)))
            
            if (size.int == size) 
              learnmatrix <- t(apply(learnmatrix, 1, 
                                     function(z) setdiff(1:n, z)))
          }
      } else {
        taby <- table(y)
        prop <- taby/sum(taby)
        siz <- n - floor(n/fold)
        classize <- .roundVector(prop * siz, siz)
        
        indlist <- sapply(names(taby), function(z) which(y == z), 
                          simplify = FALSE)
        templist <- vector(mode = "list", length = length(indlist))
        # learnmatrix <- matrix(0, niter*fold, siz) lower <- cumsum(c(1,
        # classize[-length(classize)])) upper <- cumsum(classize) 
        #templist[[i]]
        for (i in 1:length(indlist)) {
          outp <- do.call("generateLearningsets", args = list(n = taby[i], 
                                                              method = "CV", niter = niter, 
                                                              fold = fold))@learnmatrix
          templist[[i]] <- t(apply(outp, 1, function(z) ifelse(z == 0, 0,
                                                               indlist[[i]][z])))
          # learnmatrix[1:fold,lower[i]:upper[i]] <- t()) 
          #templist <- lapply(templist,
          # function(z) ifelse(z == 0, 0, indlist[[i]][z]))
        }
        
        topass <- lapply(templist, function(z) z[1:fold, , drop = FALSE])
        swaporder <- .rowSwaps(topass)
        nrep <- 1
        while (nrep < niter) {
          swaporder <- rbind(swaporder, swaporder[1:fold, , 
                                                  drop = FALSE] + fold * nrep)
          nrep <- nrep + 1
        }
        
        for (i in 1:length(templist)) templist[[i]] <- 
          templist[[i]][swaporder[, i], ]
        learnmatrix <- templist[[1]]
        for (i in 2:length(indlist)) learnmatrix <- 
          cbind(learnmatrix, templist[[i]])
      }
      } else if (method == "LOOCV") {
        learnmatrix <- matrix(rep(1:n, each = n - 1), nrow = n)
      } else if (method == "bootstrap") {
        if (is.null(niter)) 
          stop("If 'method=bootstrap', the argument 'niter' must be
               given. \n")
        
        if(!missing(y)){
          taby <- table(y)
          if (any(taby) < 1) {
            stop("There are either no censored or uncensored data. 
                 Argument 'strat' set to 'FALSE'. \n")
            strat<-FALSE	
          }}
        
        
        if (!strat) 
          learnmatrix <- t(replicate(niter, sample(n, replace = TRUE))) else {
            taby <- table(y)
            indlist <- sapply(names(taby), function(z) which(y == z), 
                              simplify = FALSE)
            learnmatrix <- matrix(nrow = niter, ncol = n)
            lower <- cumsum(c(1, taby[-length(taby)]))
            upper <- cumsum(taby)
            for (i in 1:length(indlist)) {
              learnmatrix[, lower[i]:upper[i]] <- t(replicate(niter, 
                                                              sample(indlist[[i]], taby[i], replace = TRUE)))
            }
          }
      } else {
        stop(paste("method must be one of:",
                   paste(eval(formals(generateLearningsets)$method), 
                         collapse=", "),"\n"))
      }
  
  if (strat & is.element(method, c("CV", "MCCV", "bootstrap"))) 
    method <- paste("stratified", method)
  new("LearningSets", learnmatrix = learnmatrix, method = method, 
      ntrain = ncol(learnmatrix), iter = nrow(learnmatrix))
        } 


.roundVector <- function(x, maxint){
  fx <- floor(x)
  aftercomma <- x-fx
  roundorder <- order(aftercomma, decreasing=TRUE)
  i <- 1
  while(sum(fx) < maxint){ 
    fx[roundorder[i]] <- ceiling(x[roundorder[i]])
    i <- i+1
  }
  return(fx)
}

.rowSwaps <- function(blocklist){
  cols <- length(blocklist)
  fold <- nrow(blocklist[[1]])
  learnmatrix <- blocklist[[1]]
  for(i in 2:cols) learnmatrix <- cbind(learnmatrix, blocklist[[i]])
  rs <- rowSums( learnmatrix == 0)
  allowedzeros <- ceiling(sum(rs)/fold)
  indmatrix <-  matrix(rep(1:fold, each=cols), nrow=fold, byrow=TRUE) 
  while(any(rs > allowedzeros)){
    indmatrix <- replicate(cols, sample(1:fold))
    temp2list <- blocklist
    for(i in 1:cols) temp2list[[i]] <- blocklist[[i]][indmatrix[,i], ]
    learnmatrix <- temp2list[[1]]
    for(i in 2:cols) learnmatrix <- cbind(learnmatrix, temp2list[[i]])
    rs <- rowSums( learnmatrix == 0)
  }
  return(indmatrix)
}

