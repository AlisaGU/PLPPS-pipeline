setGeneric("computeZ", function(Xlist, ylist, ...) 
  standardGeneric("computeZ"))


setMethod("computeZ", signature(Xlist = "list", ylist = "list"), 
          function(Xlist, ylist, survmethod, measure,packagetune=FALSE,add.surv = 
                     list(), add.tune = list(),add.measure=list(), 
                   plot=FALSE,trace = TRUE,names=NULL,seed=238,...) {
            set.seed(seed)
            ll <- list(...)
            if (missing(survmethod)) 
              stop("argument 'survmethod' is missing \n")
            
            if(length(Xlist)!=length(ylist))
              stop("Arguments Xlist and ylist are of different length.\n")
            
            nZ<-length(Xlist)
            
            Zmat<-matrix(NA,ncol=nZ,nrow=nZ)
            
            if(is.null(names)==TRUE)
              names=as.character(1:nZ)
            
            for(ii in 1:length(Xlist)){
              
              ###tuning
              if(packagetune==F & length(add.tune)!=0){
                
                
                add.tunet<-c(add.tune,add.surv,ll)    
                add.tunet$y<-ylist[[ii]]
                add.tunet$X<-Xlist[[ii]]
                if(is.null(add.tunet$measure)==TRUE)
                  add.tunet$measure<-measure
                add.tunet$survmethod<-survmethod
                try(tunet<-do.call(tune,add.tunet)) 
              }
              
              #fit final model
              add.survc<-c(add.surv,ll)
              add.survc$y<-ylist[[ii]]
              add.survc$X<-Xlist[[ii]]
              add.survc$survmethod<-survmethod
              if(packagetune==F & length(add.tune)!=0)
                add.survc$tuneres<-tunet
              
              
              mod<-try(do.call(learnSurvival,add.survc))
              
              if(class(mod)!='try-error'){
                
                for(jj in setdiff(1:nZ,ii)){
                  evlist<-add.measure
                  evlist$object<-mod
                  evlist$measure<-measure
                  evlist$newdata<-Xlist[[jj]]
                  evlist$newy<-ylist[[jj]]
                  val<-try(do.call(evaluate,evlist)@result[[1]])
                  if(class(val)!='try-error')
                    Zmat[ii,jj]<-val
                  
                }} else{
                  warning(paste('Error fitting model on dataset ',ii,
                                '.\n Continuing,...\n',sep=''))
                }
              
              
              
            }
            
            
            
            if(plot==TRUE){
              colnames(Zmat)<-names
              row.names(Zmat)<-names
              
              opar <- par()   
              heatmap(Zmat,scale='none',Rowv=NA,Colv=NA,col=c("#3333BB",
                                                              "#4444AA","#555599","#666688","#777777",
                                                              "#888866", "#999955","#AAAA44","#BBBB33",
                                                              "#CCCC22","#DDDD11","#EEEE00"),breaks=c(0,0.45,
                                                                                                      0.5,0.53,0.56,0.6,0.64,0.68,0.74,0.8,0.85,0.9,
                                                                                                      1),margins=c(15,15))
              
              par(fig=c(0.6,1,0,0.45),new=T)
              
              hist(unlist(Zmat),breaks=c(0,0.45,0.5,0.53,0.56,0.6,0.64,0.68,
                                         0.74,0.8,0.85,0.9,1),col=c("#3333BB","#4444AA",
                                                                    "#555599","#666688","#777777","#888866", 
                                                                    "#999955","#AAAA44","#BBBB33","#CCCC22",
                                                                    "#DDDD11","#EEEE00"),main='',xlab='Value')
              
              par(opar)
              
              
              
            }
            
            return(Zmat)
          })

