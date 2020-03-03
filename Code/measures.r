# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# fit lasso-cox with clinical covariate adjustment.
M = 5
K = c(2, 3, 5)        #2, 3, 5
fold = 4
Grid.s = 100
FIT = vector("list", M)
list.names = c("brca", "gbm", "laml", "lusc", "skcm")
names(FIT) = list.names
FIT[[1]] = FIT[[2]] = FIT[[3]] = FIT[[4]] = FIT[[5]] = vector("list", 3)

for (m in 1:M)
{
      for (p0 in 1:3)
      {
         k = K[p0]
         dat.m = DAT[[m]]
         pc = ncol(dat.m[[1]][, -(1:2)]); px = ncol(dat.m[[k]])
         penalty.factor = rep(1, pc+px); penalty.factor[1:pc] = 0
         FIT[[m]][[p0]] = cv.glmnet(as.matrix(cbind(dat.m[[1]][, -(1:2)], dat.m[[k]])), Surv(dat.m[[1]][, 1], dat.m[[1]][,2]),
                          penalty.factor = penalty.factor, family="cox", nfolds=4, nlambda=Grid.s, lambda.min.ratio=0.005)
      }
}
rm(list = ls()[-c(which(ls()=="DAT"), which(ls()=="FIT"))] )
save.image("lasso-cox.RData")





# similarity analysis

result = result1 = data.frame(matrix(0, nrow=60, ncol=11))
top = 40
COR =  vector("list", 3)
Grid.s = 100
K = c(2, 3, 5)

m = 3; n = 4


for (m in 1:5) {

 for(n in 3:5){
  if(n != m) {
   for (p0 in 1:3)
   {

      k = K[p0]


      fit.m = FIT[[m]][[p0]]
      pc.m = ncol(DAT[[m]][[1]][, -(1:2)])
      fit.n = FIT[[n]][[p0]]
      pc.n = ncol(DAT[[n]][[1]][, -(1:2)])

      # similarity at the optimal tuning point.
      # nz.m and nz.n don't include the clinic variables.
      if(sum(coef(fit.m, s=fit.m$lambda.min)[-(1:pc.m)]!=0) != 0)  nz.m = which(coef(fit.m, s=fit.m$lambda.min)[-(1:pc.m)]!=0) else
      nz.m = which(fit.m$glmnet.fit$beta[, which.min(abs(fit.m$glmnet.fit$df-top-pc.m))] != 0)[-(1:pc.m)] - pc.m
      A = data.matrix(DAT[[m]][[k]][, nz.m]); colnames(A) = colnames(DAT[[m]][[k]])[nz.m]

      if(sum(coef(fit.n, s=fit.n$lambda.min)[-(1:pc.n)]!=0) != 0)  nz.n = which(coef(fit.n, s=fit.n$lambda.min)[-(1:pc.n)]!=0) else
      nz.n = which(fit.n$glmnet.fit$beta[, which.min(abs(fit.n$glmnet.fit$df-top-pc.n))] != 0)[-(1:pc.n)] - pc.n
      B = data.matrix(DAT[[m]][[k]][, nz.n]); colnames(B) = colnames(DAT[[m]][[k]])[nz.n]


      if(length(nz.m)==0 & length(nz.n)==0) result[r, 3] = paste(0, "/", 0, sep="")   else if (length(nz.m)==0|length(nz.n)==0)
      {
         if(length(nz.m)==0) D1 = 0 else d1 = svd(A)$d; D1 = d1[d1>0.1]
         if(length(nz.n)==0) D2 = 0 else d2 = svd(B)$d; D2 = d2[d2>0.1]
         result[r, 3] = paste(0, "/", sum(svd(cbind(A,B))$d >= min(min(D1), min(D2))), sep="")
         }  else  {
         d1 = svd(A)$d; D1 = d1[d1>0.1]
         d2 = svd(B)$d; D2 = d2[d2>0.1]
         result[r, 3]  = paste((length(D1) + length(D2) - sum(svd(cbind(A,B))$d > min(min(D1), min(D2)))), "/", sum(svd(cbind(A,B))$d > min(min(D1), min(D2))), sep="")
      }

      fit0 = similar(A, B)
      result[r, 1] = length(intersect(colnames(A), colnames(B)))/length(union(colnames(A), colnames(B)))
      result[r, 2] = fit0$Rank[1];          result1[r, 3] = fit0$Rank[2]
      result[r, 4] = fit0$Cca;              result1[r, 4] = fit0$Rank[3]
      result[r, 5] = fit0$Cor[1];           result1[r, 5] = fit0$Cor[2]
      result[r ,6] = fit0$Rsqu[1];          result1[r, 6] = fit0$Rsqu[2]

 
# similarity based on solution path
      S = 100
      grid.m = grid.n = 1:S
      
      center = which.min(fit.m$cvm)
      # if the number of variables is larger than the sample size, 0 is given for that weight.
      size = min(ifelse(max(fit.m$glmnet.fit$df)<nrow(as.matrix(DAT[[m]][[k]])), length(fit.m$glmnet.fit$df), which.min(abs(fit.m$glmnet.fit$df - nrow(as.matrix(DAT[[m]][[k]]))))),
                 which.min(abs(fit.n$glmnet.fit$df - nrow(as.matrix(DAT[[n]][[k]])))))
      wei = wSum(1:S, center, size)


      Rank = Rank1 = Rank2 = Cca = Cor = Cor1 = Rsqu = Rsqu1 = Overlap =  numeric(S)
      Cors = vector("list", S)
      for(i in 1:min(which(wei==0)))
      {
        nz.m = which(fit.m$glmnet.fit$beta[, grid.m[i]]!=0)[-(1:pc.m)] - pc.m
        nz.n = which(fit.n$glmnet.fit$beta[, grid.n[i]]!=0)[-(1:pc.n)] - pc.n

        if (length(nz.m)*length(nz.n)==0)   {
           Overlap[i] = 0
           Rank[i] = 0
           Cca[i] = 0
           Cor[i] = 0
           Cors[[i]] = 0
           Rsqu[i] = 0
        } else {
           A = data.matrix(DAT[[m]][[k]][, nz.m])
           B = data.matrix(DAT[[m]][[k]][, nz.n])
           fit = similar(A, B)
           Overlap[i] = length(intersect(colnames(A), colnames(B)))/min(c(ncol(A), ncol(B)))
           Rank[i] = fit$Rank[1]; Rank1[i] = fit$Rank[2]; Rank2[i] = fit$Rank[3]
           Cca[i] = fit$Cca
           Cor[i] = fit$Cor[1]; Cor1[i] = fit$Cor[2]
           Cors[[i]] = fit$Cors
           Rsqu[i] = fit$Rsqu[1];  Rsqu1[i] = fit$Rsqu[2]
        }
        print(i)
      }

      result[r, 7] = sum((wei*Overlap))
      result[r, 8] = sum((wei*Rank));        result1[r, 8] = sum((wei*Rank1))
      result[r, 9] = sum((wei*Cca));         result1[r, 9] = sum((wei*Rank2))
      result[r, 10] = sum((wei*Cor));         result1[r, 10] = sum((wei*Cor1))
      result[r, 11] = sum((wei*Rsqu),na.rm=TRUE);       result1[r, 11] = sum((wei*Rsqu1))
      r = r + 1  
    }
   }
  }
} 

#write.csv(result[,1:6], "optimal.csv", eol = "\n")
  
