lower2symmetric<-function(n,tri)
{
  Mat<-matrix(nrow=n,ncol = n)
  Mat[lower.tri(Mat)]<-tri
  for(i in 1:(n-1))
    for(j in (i+1):n)
    {
      Mat[i,j]<-tri[(i-1)*(n-1)-(i-1)*(i-2)/2+j-i]
    }
  return(Mat)
}