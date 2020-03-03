# meta analysis --- fisher's methods.

fisher.meta = function(p) {
  chisq = -2*sum(log(p), na.rm=TRUE )
  p.value = 1-pchisq(chisq, df = 2*length(p))
  ifelse(any(p==1), 1, p.value)
}
#p <- runif(10)
#fisher.meta(p)