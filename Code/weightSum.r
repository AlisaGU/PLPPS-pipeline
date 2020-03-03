#
wSum = function(o, center, size)
{
     p =  length(o)
     w = numeric(p)
     for (i in 1:p)
     {
        if(o[i]<=center)  {
            w[i] = 2*o[i]/(center*size) 
        }  else  {
               if (o[i]>center & o[i]<=size)  
               {
                  w[i] = 2*(size-o[i])/(size*(size-center)) 
               }  else w[i] = 0
        }    
                   
     }
    # sum(w*o)
    w
}

f = function(i, center, size)
{
    
        if(i<=center)  w = 2*i/(center*size) else {
            if (i>center & i<=size) w = 2*(size-i)/(size*(size-center)) else {
            w = 0
            }
        }
             
    w
}




size = 50
o = 4*(1:25)
center= 20

jpeg("weight.jpeg")
#plot(o, wSum(o, center, size), xlab="k", ylab="", "l")
plot(1, type="n", xlim=c(0,100), ylim=c(0,0.04), xlab="k", ylab="")
lines(x=c(0,20,50,100), y=c(0,0.04,0,0))
dev.off()
ff = function(x) f(x, center, size)
integrate(ff, 0, size)

y = numeric(100)
for(i in 1:100)  y[i] = ff(i)
plot(1:100, y)



