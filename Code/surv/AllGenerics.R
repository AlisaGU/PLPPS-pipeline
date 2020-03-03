###all generics 
###(except those belonging to slot accessor or replacement functions,
### s. AllClasses.r)


setGeneric("plot")

setGeneric("predict")

setGeneric("compare", function(obj1, obj2, ...) standardGeneric("compare"))


setGeneric("evaluate", function(object, ...) standardGeneric("evaluate"))

setGeneric("aggregateMeasure", function(object, foldmeasures, ...) 
  standardGeneric("aggregateMeasure"))

setGeneric("gsaWilcoxSurv",    function(gmt,X,y,genenames,statistics,...) 
  standardGeneric("gsaWilcoxSurv"))

setGeneric("gsaTranslateGmt",  function(gmt,X,genes) 
  standardGeneric("gsaTranslateGmt"))