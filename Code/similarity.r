# measure the similarity of two matrices.
require(glmnet)
similar = function(A, B)
{
    # --------------------------------------------------------------------------
    # collinearity
    d1 = svd(A)$d; D1 = d1[d1>0.1]#[d1> max(quantile(d1, probs=0.05), 0.1)]
    d2 = svd(B)$d; D2 = d2[d2>0.1]#[d2> max(quantile(d2, probs=0.05), 0.1)]
    # Jaccard similary index
    Rank = (length(D1) + length(D2) - sum(svd(cbind(A,B))$d > min(min(D1), min(D2))))/sum(svd(cbind(A,B))$d > min(min(D1), min(D2)))
    # different tolerance values.
    D1 = d1[d1>10^-3]; D2 = d2[d2>10^-3]
    Rank1 = (length(D1) + length(D2) - sum(svd(cbind(A,B))$d > min(min(D1), min(D2))))/sum(svd(cbind(A,B))$d > min(min(D1), min(D2)))

    D1 = d1[d1>10^-6]; D2 = d2[d2>10^-6]
    Rank2 = (length(D1) + length(D2) - sum(svd(cbind(A,B))$d > min(min(D1), min(D2))))/sum(svd(cbind(A,B))$d > min(min(D1), min(D2)))


    # canical correlation analysis
    Cca = sum(cancor(A,B)$cor > 0.5)/min(c(ncol(A), ncol(B)))

    # correlation values which are seem to be outlier.
    rho = cor(A, B)
    acor = abs(rho)
    cutoff = (exp(2*1.96/sqrt(nrow(A)-3))-1)/(exp(2*1.96/sqrt(nrow(B)-3))+1)
    Cor = sum(acor>cutoff)/prod(c(ncol(A), ncol(B)))
    # just use the mean.
    Cor1 = mean(acor)

    p.m = ncol(A); n = nrow(A); index = ceiling(sample(1:n)/n * 2) 
    Rsqu.A = numeric(p.m)
    
    for (v in 1:p.m)
    {
        if (ncol(B)==1)   {
           fit = lm(A[index!=1, v]~B[index!=1, ])  
           Rsqu.A[v] = summary(fit)$r.squared
        } else {
           fit = cv.glmnet(B[index!=1, ], A[index!=1, v])  
           nz = which(coef(fit, s="lambda.min")[-1]!=0)
           if(length(nz)==0)  Rsqu.A[v] = 0 else Rsqu.A[v] = summary(lm(A[index==1,v]~ B[index==1, nz]))$r.squared
        }   
    }
      
    Rsqu = sum(Rsqu.A>0.5)/p.m
    Rsqu1 = mean(Rsqu.A)

    list(Rank=c(Rank, Rank1, Rank2), Cca=Cca, Cor=c(Cor, Cor1), Cors=c(rho), Rsqu=c(Rsqu, Rsqu1))
}
