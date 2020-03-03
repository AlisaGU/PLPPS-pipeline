setwd("/TCGA_Similarity")
#tumor: "brca", "gbm", "laml", "lusc", "skcm"
# clinic, mRNA, methylation, miRNA, CNA
source("impute.r")
source("overlap.r")
source("fisher.meta.r")
source("pvalue.r")
source("similarity.r")
require(glmnet)
require(survival)

filenames = list.files("/TCGA_Similarity/data")
setwd("/TCGA_Similarity/data")
ldf = lapply(filenames, read.table, header=TRUE)
set.seed(2013201414)
setwd("/TCGA_Similarity")


# brca
brca = list(ldf[[1]], ldf[[5]], ldf[[3]], ldf[[4]], ldf[[2]])
brca.samples = Reduce(intersect, list(rownames(brca[[1]][-union(which(brca[[1]][,1]==0), which(is.na(brca[[1]][,1])) ) , ]),
                                      rownames(brca[[2]]),
                                      rownames(brca[[3]]),
                                      rownames(brca[[4]]),
                                      rownames(brca[[5]])))
                                      
Brca = list(brca[[1]][match(brca.samples, rownames(brca[[1]])), ],
            impute(brca[[2]][match(brca.samples, rownames(brca[[2]])), ]),
            impute(brca[[3]][match(brca.samples, rownames(brca[[3]])), ]),
            brca[[4]][match(brca.samples, rownames(brca[[4]])), ],
            brca[[5]][match(brca.samples, rownames(brca[[5]])), ])

                                      
                                      
   
# gbm, no missing in the survivla times
gbm = list(ldf[[6]], ldf[[10]], ldf[[8]], ldf[[9]], ldf[[7]])
gbm.samples = Reduce(intersect, list(rownames(gbm[[1]]),
                                      rownames(gbm[[2]][-which(rowSums(is.na(gbm[[2]]))==5178), ]),
                                      rownames(gbm[[3]]),
                                      rownames(gbm[[4]]),
                                      rownames(gbm[[5]])))
                                      
Gbm = list(gbm[[1]][match(gbm.samples, rownames(gbm[[1]])), ],
           impute(gbm[[2]][match(gbm.samples, rownames(gbm[[2]])), ]),
           impute(gbm[[3]][match(gbm.samples, rownames(gbm[[3]])), ]),
           gbm[[4]][match(gbm.samples, rownames(gbm[[4]])), ],
           gbm[[5]][match(gbm.samples, rownames(gbm[[5]])), ])


# laml
laml = list(ldf[[11]], ldf[[14]], ldf[[13]], NULL, ldf[[12]])
laml.samples = Reduce(intersect, list(rownames(laml[[1]][-union(which(laml[[1]][,1]==0), which(is.na(laml[[1]][,1])) ) , ]),
                                      rownames(laml[[2]]),
                                      rownames(laml[[3]]),
                                      #rownames(laml[[4]]),
                                      rownames(laml[[5]])))

Laml = list(laml[[1]][match(laml.samples, rownames(laml[[1]])), ],
            laml[[2]][match(laml.samples, rownames(laml[[2]])), ],
            impute(laml[[3]][match(laml.samples, rownames(laml[[3]])), ]),
            NULL,
            laml[[5]][match(laml.samples, rownames(laml[[5]])), ])


# lusc
lusc = list(ldf[[15]], ldf[[19]], ldf[[17]], ldf[[18]], ldf[[16]])
lusc.samples = Reduce(intersect, list(rownames(lusc[[1]][-union(which(lusc[[1]][,1]==0), which(is.na(lusc[[1]][,1])) ) , ]),
                                      rownames(lusc[[2]]),
                                      rownames(lusc[[3]]),
                                      rownames(lusc[[4]]),
                                      rownames(lusc[[5]])))

Lusc = list(lusc[[1]][match(lusc.samples, rownames(lusc[[1]])), ],
            impute(lusc[[2]][match(lusc.samples, rownames(lusc[[2]])), ]),
            impute(lusc[[3]][match(lusc.samples, rownames(lusc[[3]])), ]),
            lusc[[4]][match(lusc.samples, rownames(lusc[[4]])), ],
            lusc[[5]][match(lusc.samples, rownames(lusc[[5]])), ])


# melanoma
mela = list(ldf[[20]][,-(6:7)], ldf[[23]], ldf[[22]], NULL, ldf[[21]])
# sample's id in mela[[2]] are different with others.
rownames(mela[[2]]) =  paste(substr(rownames(mela[[2]]), 1, 4), substr(rownames(mela[[2]]), 6, 7), substr(rownames(mela[[2]]), 9, 12), sep="-")
mela.samples = Reduce(intersect, list(rownames(mela[[1]][-union(which(mela[[1]][,1]==0), which(is.na(mela[[1]][,1])) ) , ]),
                                     rownames(mela[[2]]),
                                      rownames(mela[[3]]),
                                      #rownames(mela[[4]]),                                           
                                      rownames(mela[[5]])))



Mela = list(mela[[1]][match(mela.samples, rownames(mela[[1]])), ],
            mela[[2]][match(mela.samples, rownames(mela[[2]])), ],
            impute(mela[[3]][match(mela.samples, rownames(mela[[3]])), ]),
            NULL,
            # remove markers with lots of missings.
            impute(mela[[5]][match(mela.samples, rownames(mela[[5]])), ][, -which(colSums(is.na(mela[[5]][match(mela.samples, rownames(mela[[3]])), ]))==338)]))

   
# ------------------------------------------------------------------------------
mrna.markers = Reduce(intersect, list(colnames(Brca[[2]]),
                                      colnames(Gbm[[2]]),
                                      colnames(Laml[[2]]),
                                      colnames(Lusc[[2]]),
                                      colnames(Mela[[2]])))
                                      
methylation.markers = Reduce(intersect, list(colnames(Brca[[3]]),
                                             colnames(Gbm[[3]]),
                                             colnames(Laml[[3]]),
                                             colnames(Lusc[[3]]),
                                             colnames(Mela[[3]])))
                                              
mirna.markers = Reduce(intersect, list(colnames(Brca[[4]]),
                                       colnames(Gbm[[4]]),
                                       colnames(Laml[[4]]),
                                       colnames(Lusc[[4]]),
                                       colnames(Mela[[4]])))

cna.markers = Reduce(intersect, list(colnames(Brca[[5]]),
                                     colnames(Gbm[[5]]),
                                     colnames(Laml[[5]]),
                                     colnames(Lusc[[5]]),
                                     colnames(Mela[[5]])))

BRCA = list(Brca[[1]], Brca[[2]][, match(mrna.markers, colnames(Brca[[2]]))], 
                       Brca[[3]][, match(methylation.markers, colnames(Brca[[3]]))], 
                       Brca[[4]][, match(mirna.markers, colnames(Brca[[4]]))], 
                       Brca[[5]][, match(cna.markers, colnames(Brca[[5]]))])
                       
GBM = list(Gbm[[1]], Gbm[[2]][, match(mrna.markers, colnames(Gbm[[2]]))], 
                       Gbm[[3]][, match(methylation.markers, colnames(Gbm[[3]]))], 
                       Gbm[[4]][, match(mirna.markers, colnames(Gbm[[4]]))], 
                       Gbm[[5]][, match(cna.markers, colnames(Gbm[[5]]))])
                         
LAML = list(Laml[[1]], Laml[[2]][, match(mrna.markers, colnames(Laml[[2]]))], 
                       Laml[[3]][, match(methylation.markers, colnames(Laml[[3]]))], 
                       Laml[[4]][, match(mirna.markers, colnames(Laml[[4]]))], 
                       Laml[[5]][, match(cna.markers, colnames(Laml[[5]]))])

LUSC = list(Lusc[[1]], Lusc[[2]][, match(mrna.markers, colnames(Lusc[[2]]))], 
                       Lusc[[3]][, match(methylation.markers, colnames(Lusc[[3]]))], 
                       Lusc[[4]][, match(mirna.markers, colnames(Lusc[[4]]))], 
                       Lusc[[5]][, match(cna.markers, colnames(Lusc[[5]]))])
                       
MELA = list(Mela[[1]], Mela[[2]][, match(mrna.markers, colnames(Mela[[2]]))], 
                       Mela[[3]][, match(methylation.markers, colnames(Mela[[3]]))], 
                       Mela[[4]][, match(mirna.markers, colnames(Mela[[4]]))], 
                       Mela[[5]][, match(cna.markers, colnames(Mela[[5]]))])
                       


sum(is.na(BRCA[[1]])); sum(is.na(BRCA[[2]])); sum(is.na(BRCA[[3]])); sum(is.na(BRCA[[5]])); 
sum(is.na(GBM[[1]])); sum(is.na(GBM[[2]])); sum(is.na(GBM[[3]])); sum(is.na(GBM[[5]])); 
sum(is.na(LAML[[1]])); sum(is.na(LAML[[2]])); sum(is.na(LAML[[3]])); sum(is.na(LAML[[5]])); 
sum(is.na(LUSC[[1]])); sum(is.na(LUSC[[2]])); sum(is.na(LUSC[[3]])); sum(is.na(LUSC[[5]])); 
sum(is.na(MELA[[1]])); sum(is.na(MELA[[2]])); sum(is.na(MELA[[3]])); sum(is.na(MELA[[5]])); 

# overlaps of methylations are too small.
dim(BRCA[[1]]); dim(BRCA[[2]]); dim(BRCA[[3]]); dim(BRCA[[4]]); dim(BRCA[[5]])   
dim(GBM[[1]]); dim(GBM[[2]]); dim(GBM[[3]]); dim(GBM[[4]])  ; dim(GBM[[5]])                                               
dim(LAML[[1]]); dim(LAML[[2]]); dim(LAML[[3]]); dim(LAML[[4]]); dim(LAML[[5]])                                                                                           
dim(LUSC[[1]]); dim(LUSC[[2]]); dim(LUSC[[3]]); dim(LUSC[[4]]); dim(LUSC[[5]])    
dim(MELA[[1]]); dim(MELA[[2]]); dim(MELA[[3]]); dim(MELA[[4]]); dim(MELA[[5]])    


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# focus on mRNA, methylation and CNA first, screening variables to get 2500 markers with meta-pvalues.
Dat = list(BRCA, GBM, LAML, LUSC, MELA)
M = 5
p = 2500
#nz = vector(mode="list", length=M)

Index = vector(mode="list", length=5)
k = 2 # if k==3, adjust index.
for (k in c(2, 3, 5))
{
  if (k==3) Index[[k]] = 1:ncol(Dat[[1]][[k]]) else {
      P = ncol(Dat[[1]][[k]]); Pvalues = matrix(0, ncol=M, nrow=P)
      for (m in 1:M)
          Pvalues[, m] = pvalue(Dat[[m]][[k]], Surv(Dat[[m]][[1]][, 1], Dat[[m]][[1]][,2]))
          pool = numeric(P)
      for (i in 1:P)
          pool[i] = fisher.meta(Pvalues[i, ]) 
      Index[[k]] = order(pool, decreasing = FALSE)[1:p]
  }
}

DAT = list(    # clinic,    mRNA,                     methylation,          miRNA,        CNA           
   BRCA = list(BRCA[[1]], BRCA[[2]][, Index[[2]]], BRCA[[3]][, Index[[3]]], NULL, BRCA[[5]][, Index[[5]]]),
   GBM = list(GBM[[1]], GBM[[2]][, Index[[2]]], GBM[[3]][, Index[[3]]], NULL, GBM[[5]][, Index[[5]]]),
   LAML = list(LAML[[1]], LAML[[2]][, Index[[2]]], LAML[[3]][, Index[[3]]], NULL, LAML[[5]][, Index[[5]]]),
   LUSC = list(LUSC[[1]], LUSC[[2]][, Index[[2]]], LUSC[[3]][, Index[[3]]], NULL, LUSC[[5]][, Index[[5]]]),
   SKCM = list(MELA[[1]], MELA[[2]][, Index[[2]]], MELA[[3]][, Index[[3]]], NULL, MELA[[5]][, Index[[5]]]))

for (i in 1:5)
    for (j in 1:5)
          print(dim(DAT[[i]][[j]]))
rm(list = ls()[-which(ls()=="DAT")] )
#save.image("fiveCancers.RData")



top = 40
# without covariates  adjustment.
for (m in 1:M)
{
   dat = DAT[[m]]
   #                         predictors                         y,            delta   
   fit = cv.glmnet(as.matrix(dat[[k]][, index]), Surv(dat[[1]][, 1], dat[[1]][,2]), family="cox", nfolds=4)
   nz[[m]] = which(fit$glmnet.fit$beta[, which.min(abs(fit$glmnet.fit$df-top))] != 0) 
}
overlap(nz)

# with covariates  adjustment.          
for (m in 1:M)
{
   dat = Dat[[m]]
   penalty.factor = rep(1, ncol(dat[[k]][, index])+ncol(dat[[1]][,-(1:2)]))
   penalty.factor[1:ncol(dat[[1]][,-(1:2)])] = 0
   #                         predictors                         y,            delta   
   fit = cv.glmnet(as.matrix(cbind(dat[[1]][,-(1:2)] , dat[[k]][, index])), Surv(dat[[1]][, 1], dat[[1]][,2]), family="cox", penalty.factor=penalty.factor)
   nz[[m]] = which(fit$glmnet.fit$beta[-(1:ncol(dat[[1]][,-(1:2)])), which.min(abs(fit$glmnet.fit$df-ncol(dat[[1]][,-(1:2)])-top))] != 0) 
}
overlap(nz)

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# prediction results without clinical covariate adjustment .
M = 5
K = c(2, 3, 5)        #2, 3, 5
fold = 4
top = 40
sim.n = 50
result = result.sd = array(0, dim=c(M, M, 3))
for (p0 in 1:3)
{
   k = K[p0]
   logrank = array(0, dim=c(M, M, sim.n))
   for (sim in 1:sim.n)
   {
      for (m in 1:M)
      { 
         dat.m = DAT[[m]]
         index.m = ceiling(sample(1:nrow(dat.m[[1]]))/nrow(dat.m[[1]]) * fold)       
         fit = cv.glmnet(as.matrix(dat.m[[k]][which(index.m!=1), ]), Surv(dat.m[[1]][which(index.m!=1), 1], dat.m[[1]][which(index.m!=1),2]), family="cox", nfolds=4)
      
         # prediction based on markers selected by itself.
         nz.m = which(fit$glmnet.fit$beta[, which.min(abs(fit$glmnet.fit$df-top))] != 0) 
         fit.m = coxph(Surv(dat.m[[1]][which(index.m!=1), 1], dat.m[[1]][which(index.m!=1),2]) ~ ridge(as.matrix(dat.m[[k]][which(index.m!=1), ])[, nz.m],  theta=10, scale=TRUE) ) 
         na = which(is.na(coef(fit.m)))
         if(length(na)!=0) lp.m = c(as.matrix(dat.m[[k]][which(index.n==1), ][, nz.m][,-na]) %*% coef(fit.m)[-na])  else
         lp.m = c(as.matrix(dat.m[[k]][which(index.m==1), ][, nz.m]) %*% coef(fit.m))
         
         group.m = numeric(sum(index.m==1)); group.m[which(lp.m >= median(lp.m))] = 1 
         logrank[m,m,sim] = ifelse(sd(group.m)==0, 0, survdiff(Surv(dat.m[[1]][which(index.m==1), 1], dat.m[[1]][which(index.m==1),2]) ~ group.m)$chisq)
   
         # prediction based on markers selected by others. 
         for (n in 1:M)
         {  
            if(n!=m)
            { 
               dat.n = DAT[[n]]
               index.n = ceiling(sample(1:nrow(dat.n[[1]]))/nrow(dat.n[[1]]) * fold)    
               fit.n = coxph(Surv(dat.n[[1]][which(index.n!=1), 1], dat.n[[1]][which(index.n!=1),2]) ~ ridge(as.matrix(dat.n[[k]][which(index.n!=1), ])[, nz.m],  theta=10, scale=TRUE) ) 
               na = which(is.na(coef(fit.n)))
               # some coefficient go to infinity.
               if(length(na)!=0) lp.n = c(as.matrix(dat.n[[k]][which(index.n==1), ][, nz.m][,-na]) %*% coef(fit.n)[-na])  else
               lp.n = c(as.matrix(dat.n[[k]][which(index.n==1), ][, nz.m]) %*% coef(fit.n))
            
               group.n = numeric(sum(index.n==1)); group.n[which(lp.n >= median(lp.n))] = 1
               logrank[m,n,sim] = ifelse(sd(group.n)==0, 0, survdiff(Surv(dat.n[[1]][which(index.n==1), 1], dat.n[[1]][which(index.n==1),2]) ~ group.n)$chisq)
             } 
          }  
      }
   }
   
   result[, , p0] = apply(logrank, c(1,2), mean)
   result.sd[,,p0] = apply(logrank, c(1,2), sd)
}


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# prediction results with clinical covariate adjustment and these variables are used for prediction.
M = 5
K = c(2, 3, 5)        #2, 3, 5
fold = 4
top = 40
sim.n = 50
result = result.sd = array(0, dim=c(M, M, 3))
for (p0 in 1:3)
{
   k = K[p0]
   logrank = array(0, dim=c(M, M, sim.n))
   for (sim in 36:sim.n)
   {
      for (m in 1:M)
      { 
         dat.m = DAT[[m]]
         index.m = ceiling(sample(1:nrow(dat.m[[1]]))/nrow(dat.m[[1]]) * fold)  
         pc = ncol(dat.m[[1]][, -(1:2)]); px = ncol(dat.m[[k]])       
         penalty.factor = rep(1, pc+px); penalty.factor[1:pc] = 0
         fit = cv.glmnet(as.matrix(cbind(dat.m[[1]][which(index.m!=1), -(1:2)], dat.m[[k]][which(index.m!=1), ])), Surv(dat.m[[1]][which(index.m!=1), 1], dat.m[[1]][which(index.m!=1),2]), 
                         penalty.factor = penalty.factor, family="cox", nfolds=4)
      
         # prediction based on markers selected by itself.
         # selected markers or top 40 markers, the nz.m is all variables with nonzero coefficients, including the clinic variables.
         if(sum(coef(fit, s="lambda.min")[-(1:pc)]!=0)!=0) nz.m = which(coef(fit, s="lambda.min")!=0) else
         nz.m = which(fit$glmnet.fit$beta[, which.min(abs(fit$glmnet.fit$df-top-pc))] != 0) 
         fit.m = coxph(Surv(dat.m[[1]][which(index.m!=1), 1], dat.m[[1]][which(index.m!=1),2]) ~ 
                       ridge(as.matrix(cbind(dat.m[[1]][which(index.m!=1), -(1:2)], dat.m[[k]][which(index.m!=1), ]))[, nz.m],  theta=10, scale=TRUE) ) 
         na = which(is.na(coef(fit.m)))
         if(length(na)!=0) lp.m = c(as.matrix(cbind(dat.m[[1]][which(index.m==1), -(1:2)], dat.m[[k]][which(index.m==1), ]))[, nz.m][, -na] %*% coef(fit.m)[-na])   else
         lp.m = c(as.matrix(cbind(dat.m[[1]][which(index.m==1), -(1:2)], dat.m[[k]][which(index.m==1), ]))[, nz.m] %*% coef(fit.m))
         
         
         group.m = numeric(sum(index.m==1)); group.m[which(lp.m >= median(lp.m))] = 1 
         logrank[m,m,sim] = ifelse(sd(group.m)==0, 0, survdiff(Surv(dat.m[[1]][which(index.m==1), 1], dat.m[[1]][which(index.m==1),2]) ~ group.m)$chisq)
   
         # prediction based on markers selected by others. 
         for (n in 1:M)
         {  
            if(n!=m)
            { 
               dat.n = DAT[[n]]
               index.n = ceiling(sample(1:nrow(dat.n[[1]]))/nrow(dat.n[[1]]) * fold)                    
               fit.n = coxph(Surv(dat.n[[1]][which(index.n!=1), 1], dat.n[[1]][which(index.n!=1),2]) ~ 
                             ridge(as.matrix(cbind(dat.n[[1]][which(index.n!=1), -(1:2)], dat.n[[k]][which(index.n!=1), nz.m[-(1:pc)]-pc])),  theta=10, scale=TRUE) ) 
               na = which(is.na(coef(fit.n)))
               # some coefficient go to infinity.
               if(length(na)!=0) lp.n = c(as.matrix(cbind(dat.n[[1]][which(index.n==1), -(1:2)], dat.n[[k]][which(index.n==1), nz.m[-(1:pc)]-pc])[,-na]) %*% coef(fit.n)[-na])  else
               lp.n = c(as.matrix(cbind(dat.n[[1]][which(index.n==1), -(1:2)], dat.n[[k]][which(index.n==1), nz.m[-(1:pc)]-pc])) %*% coef(fit.n))
            
               group.n = numeric(sum(index.n==1)); group.n[which(lp.n >= median(lp.n))] = 1
               logrank[m,n,sim] = ifelse(sd(group.n)==0, 0, survdiff(Surv(dat.n[[1]][which(index.n==1), 1], dat.n[[1]][which(index.n==1),2]) ~ group.n)$chisq)
             } 
          }  
      }
   }
   
   result[, , p0] = apply(logrank, c(1,2), mean)
   result.sd[,,p0] = apply(logrank, c(1,2), sd)
}


