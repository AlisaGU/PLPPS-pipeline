### signature LearnOut,LearnOut
setMethod("compare", signature(obj1 = "LearnOut", obj2 = "LearnOut"), 
          function(obj1, obj2, measure, newdata = NULL, newy = NULL, ...) {
            ### check whether LearningSets are identical
            if (any(obj1@LearningSets@learnmatrix == 
                    obj2@LearningSets@learnmatrix) == FALSE) 
              stop("Learningsets are not identical")
            complist <- list()
            
            ls <- obj1@LearningSets@learnmatrix
            for (i in 1:nrow(obj1@LearningSets@learnmatrix)) {
              if (is.null(newdata) == TRUE & is.null(newy) == TRUE) 
                complist[[i]] <- compare(obj1@ModelLearnedlist[[i]], 
                                         obj2@ModelLearnedlist[[i]], 
                                         newdata = obj1@X[-ls[i, ], ], 
                                         newy = obj1@y[-ls[i, ]], measure=measure, ...)
              else if (is.null(newdata) == FALSE & is.null(newy) == FALSE) 
                complist[[i]] <- compare(obj1@ModelLearnedlist[[i]], 
                                         obj2@ModelLearnedlist[[i]], newdata = newdata, 
                                         newy = newy, measure=measure, ...) 
              else stop("Either both newdata and newy have to be provided
                        or none of these arguments.")
            }
            
            return(complist)
            })

### signature ModelLearned,ModelLearned
setMethod("compare", signature(obj1 = "ModelLearned", obj2 = "ModelLearned"), 
          function(obj1, obj2, measure, newdata, newy, ...) {
            
            
            lp1 <- predict(obj1, newdata = newdata, type = "lp")
            lp2 <- predict(obj2, newdata = newdata, type = "lp")
            
            compare(obj1=lp1,obj2=lp2,measure=measure,newy=newy,...)	
            
          })

### signature ModelBase, ModelBase
setMethod("compare", signature(obj1 = "ModelBase", obj2 = "ModelBase"), 
          function(obj1, obj2, measure, newdata, newy, ...) {
            
            
            lp1 <- predict(obj1, newdata = newdata, type = "lp")
            lp2 <- predict(obj2, newdata = newdata, type = "lp")
            
            compare(obj1=lp1,obj2=lp2,measure=measure,newy=newy,...)	
            
          })			



###signature LinearPrediction, LinearPrediction			
setMethod("compare", signature(obj1 = "LinearPrediction", 
                               obj2 = "LinearPrediction"), function(obj1, obj2, measure, 
                                                                    newy, ...) {			
                                 
                                 obj<-new(measure)
                                 compare(obj1=obj,lp1=obj1,lp2=obj2,newy=newy,...)
                               })




###lowest level compare functions
###signature IDI
setMethod("compare",signature(obj1="IDI",obj2="missing"),
          function(obj1,obj2,lp1,lp2,newy,...){
            require(survIDINRI)
            ll<-list(...)
            ll$indata <- newy
            ll$covs0<-lp1@lp
            ll$covs1<-lp2@lp
            
            out <- do.call("IDI.INF", args = ll)
            new("IDI",estimate=out$m1[1],conf.int=out$m1[2:3])
          })

###signature NRI
setMethod("compare",signature(obj1="NRI",obj2="missing"),
          function(obj1,obj2,lp1,lp2,newy,...){
            require(survIDINRI)
            ll<-list(...)
            ll$indata <- newy
            ll$covs0<-lp1@lp
            ll$covs1<-lp2@lp
            
            out <- do.call("IDI.INF", args = ll)
            new("NRI",estimate=out$m2[1],conf.int=out$m2[2:3])
          })

###signature MIRS
setMethod("compare",signature(obj1="MIRS",obj2="missing"),
          function(obj1,obj2,lp1,lp2,newy,...){
            require(survIDINRI)
            ll<-list(...)
            ll$indata <- newy
            ll$covs0<-lp1@lp
            ll$covs1<-lp2@lp
            
            out <- do.call("IDI.INF", args = ll)
            new("MIRS",estimate=out$m3[1],conf.int=out$m3[2:3])
          })

####signature CDelta
setMethod("compare",signature(obj1="CDelta",obj2="missing"),
          function(obj1,obj2,lp1,lp2,newy,...){
            require(survC1)
            ll<-list(...)
            ll$mydata <- newy[,1:2]
            ll$covs0<-lp1@lp
            ll$covs1<-lp2@lp
            
            out <- do.call("Inf.Cval.Delta", args = ll)
            new("CDelta",estimate=out[3,1],conf.int=out[3,3:4])
          })
