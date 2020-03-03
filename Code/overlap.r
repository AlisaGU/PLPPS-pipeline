# overlap is a function used to calculate the length of intersetion of each row
# in the matrix A.

overlap = function(A)
{
     count = function(x, y, A)
     {
         if (is.matrix(A))  
            length(intersect(A[x,], A[y,])) 
         else if (is.list(A)) 
            length(intersect(A[[x]], A[[y]]))
         else print("A need to be a matrix or a list")   
     }
     m = ifelse(is.matrix(A),nrow(A),length(A))
     result = matrix(0, nrow=m, ncol=m)
     if  (is.matrix(A))  diag(result) = ncol(A) else 
         diag(result) =  sapply(A,length)
     for (i in 1:(m-1))
     {                                 
         for (j in (i+1):m)
         {
             result[i, j] = result[j, i] = count(i, j, A)
         }
     }
     result
}