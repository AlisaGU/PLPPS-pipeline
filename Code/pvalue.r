# screening variable based on cox model
# input:
#         y has form:  Surv(y, status)
#         x
require(survival)
pvalue = function(x, y)
{

    x = data.matrix(x)
    p = ncol(x)
    pvalues = numeric(p)
    
    for (i in 1:p)
    {
        if (sum(abs(x[, i]), na.rm=T)==0 | sd(x[, i], na.rm=T)==0)
        {
           pvalues[i] = 1
        }  else {
           fit = coxph(y ~ x[,i])
           pvalues[i] = summary(fit)$coefficients[1,5]
        }
    }
    pvalues

}
