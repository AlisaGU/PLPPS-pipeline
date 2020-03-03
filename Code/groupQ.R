groupQ<-function(x,Q1score,medianScore,Q3score){
  Group<-""
  if(x<=Q1score)
    Group<-"Class one"
  if(x>=Q1score&x<=medianScore)
    Group<-"Class two"
  if(x>=medianScore&x<=Q3score)
    Group<-"Class three"
  if(x>=Q3score)
    Group<-"Class four"
  return(Group)
}