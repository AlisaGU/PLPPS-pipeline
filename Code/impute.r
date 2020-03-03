# the function used to impute the missing values by the median of corresponding columns.
# a normal-like dataframe
#       x1, x2, x3 ...
# obs1   1   2   3
# obs2   2   3  NA
# ...   NA  NA   3
# 7/11/2013 in yale
# Example:
#         A = data.frame(x1=c(1,2,NA),x2=c(2,3,NA),x3=c(3,NA,3))
#         impute(A)
impute = function(dat)
{
    index = which(is.na(dat), arr.ind=T)
    
    if (length(index)!=0)
    {
       # why use loop? some times the data sizes are too large to use apply
       for (r in 1:nrow(index))
       {
            dat[index[r,1], index[r,2]] = median(dat[,index[r,2]], na.rm=TRUE)
       }
    }
    
    dat
}
