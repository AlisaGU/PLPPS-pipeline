z.test<-function(x,mu,sigma,alternative="two.sided"){
  n=length(x)
  result=list()
  mean=mean(x)
  z=(mean-mu)/(sigma/sqrt(n))
  #options(digits=4) 
  result$mean=mean;result$z=z 
  result$P=2*pnorm(abs(z),lower.tail=FALSE) 
  if(alternative=="greater") result$P=pnorm(fz,lower.tail=FALSE)
  else if(alternative=="less") result$P=pnorm(z)
  return(result)
}
