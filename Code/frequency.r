library(glmnet)
library(survival)
library(ggplot2)

# Define function 
# Function1-1 gen.freq: to generate each feature's frequency
##For measurements mRNA and CNA, the number of features after screen is 2500, using function gen.freq()
gen.freq <- function(A){ #A is glmnet object
  obj <- A$glmnet.fit
  Beta <- obj$beta
  #omit the clinical data. only fetch the 2500 features
  freq <- as.data.frame(t(as.matrix(abs(sign(Beta)))))
  freq <- freq[,((length(freq)-2499):length(freq))] 
  freq <- rbind(freq,apply(freq, 2, sum))
  freq <- sort(freq[nrow(freq),],decreasing=T)
  freq <- freq[,freq>0]
  return(freq)
}

# Function1-2 gen.freq2: to generate each feature's frequency
##For measurement methylation, the number of features after screen is 193, using function gen.freq2()
gen.freq2 <- function(A){
  obj <- A$glmnet.fit
  Beta <- obj$beta
  #omit the clinical data. only fetch the 193 features
  freq <- as.data.frame(t(as.matrix(abs(sign(Beta)))))
  freq <- freq[,(length(freq)-192):length(freq)]
  freq <- rbind(freq,apply(freq, 2, sum))
  freq <- sort(freq[nrow(freq),],decreasing=T)
  freq <- freq[,freq>0]
  return(freq)
}


# Function2-1 Overlap: calculate overlap score, corresponding to gen.freq
Overlap <- function(A, B){
  C <- list()
  freq1 <- gen.freq(A)
  freq2 <- gen.freq(B)
  print(c(length(freq1),length(freq2)))
  
  list1 <- names(freq1)
  list2 <- names(freq2)
  w1 <- match(list1, list2)
  
  C$freq1 <- freq1
  C$freq2 <- freq2
  C$raw.overlaps <- sum(abs(sign(w1)), na.rm=T)
  print(c("raw.overlaps=", C$raw.overlaps))
  
  a <- C$raw.overlaps
  total <- nrow(as.data.frame(t(merge(freq2,freq1,all=TRUE))))
  C$jaccard <- a/total
  print(c("Jaccard index=", C$jaccard))
  
  w1 <- 1- ((w1-1)/length(list2))
  w2 <- 1:length(list1)
  w2 <- 1- ((w2-1)/length(list1))
  C$overlaps <- sum(w1*w2, na.rm=T)
  print(c("overlaps=", C$overlaps))
  return(C)
}

# Function2-2 Overlap2: calculate overlap score, corresponding to gen.freq2
Overlap2 <- function(A, B){
  C <- list()
  freq1 <- gen.freq2(A)
  freq2 <- gen.freq2(B)
  print(c(length(freq1),length(freq2)))
  
  list1 <- names(freq1)
  list2 <- names(freq2)
  w1 <- match(list1, list2)
  
  C$freq1 <- freq1
  C$freq2 <- freq2
  C$raw.overlaps <- sum(abs(sign(w1)), na.rm=T)
  print(c("raw.overlaps=", C$raw.overlaps))
  
  a <- C$raw.overlaps
  total <- nrow(as.data.frame(t(merge(freq2,freq1,all=TRUE))))
  C$jaccard <- a/total
  print(c("Jaccard index=", C$jaccard))
  
  w1 <- 1- ((w1-1)/length(list2))
  w2 <- 1:length(list1)
  w2 <- 1- ((w2-1)/length(list1))
  C$overlaps <- sum(w1*w2, na.rm=T)
  print(c("overlaps=", C$overlaps))
  return(C)
}


#Data analysis processing
#Generate raw overlaps, Jaccard index and overlap score
#j denotes the first cancer type and k denotes the second one. The results are showed as an example below, 
#which means for the measure mRNA of BRCA and GBM, the numbers of features selected in 500 times are 102 
#and 285, respectively. The number of raw overlaps of them is 12, and Jaccard index is 0.032. The final 
#overlap score is about 3.097.
#[1] 1
#[1] 2
#[1] 102 285
#[1] "raw.overlaps=" "12"           
#[1] 102 285
#[1] "raw.overlaps=" "12"           
#[1] "Jaccard index=" "0.032"         
#[1] "overlaps="        "3.09680082559339"


##mRNA
for (j in 1:4){
   for (k in (j+1):5){
     print(j); print(k)
     Overlap(FIT[[j]][[1]],FIT[[k]][[1]]) #using Overlap()
   }
}

#methylation
for (j in 1:4){
  for (k in (j+1):5){
    print(j); print(k)
    Overlap2(FIT[[j]][[2]],FIT[[k]][[2]]) #using Overlap2()
  }
}

#CNA
for (j in 1:4){
  for (k in (j+1):5){
    print(j); print(k)
    Overlap(FIT[[j]][[3]],FIT[[k]][[3]]) #using Overlap()
  }
}

# Figure A3
# FunctionA3-1 Draw.overlap: generate figures of freqency, for mRNA and CNA
Draw.overlap <- function(i,j,k){
  C <- list()
  overlap.obj <- Overlap(FIT[[i]][[k]],FIT[[j]][[k]])
  freq1 <- (overlap.obj$freq1)
  freq1 <- freq1/max(freq1)
  freq2 <- (overlap.obj$freq2)
  freq2 <- freq2/max(freq2)
  
  Freq <- as.data.frame(t(merge(freq2,freq1,all=TRUE)))
  names(Freq) <- c("BRCA","GBM")
  Freq$BRCA[which(is.na(Freq$BRCA))] <- 0
  Freq$GBM[which(is.na(Freq$GBM))] <- 0
  Freq$similar <- "Low"
  Freq$similar[which(Freq$BRCA > 0.5 & Freq$GBM >0.5)] <- "High"
  print(summary(as.factor(Freq$similar)))
  C$Freq <- Freq
  high <- Freq[which(Freq$BRCA > 0.5 & Freq$GBM >0.5),]
  
  C$high <- high
  C$plot <- (qplot(BRCA,GBM, data=Freq, color=similar) +
               geom_text(aes(label=row.names(high)), position="jitter",color="black",data=high))  
  return(C)
}

# FunctionA3-2 Draw.overlap: generate figures of freqency, for methylation
Draw.overlap2 <- function(i,j,k){
  C <- list()
  overlap.obj <- Overlap2(FIT[[i]][[k]],FIT[[j]][[k]])
  freq1 <- (overlap.obj$freq1)
  freq1 <- freq1/max(freq1)
  freq2 <- (overlap.obj$freq2)
  freq2 <- freq2/max(freq2)
  
  Freq <- as.data.frame(t(merge(freq2,freq1,all=TRUE)))
  names(Freq) <- c("BRCA","GBM")
  Freq$BRCA[which(is.na(Freq$BRCA))] <- 0
  Freq$GBM[which(is.na(Freq$GBM))] <- 0
  Freq$similar <- "Low"
  Freq$similar[which(Freq$BRCA > 0.5 & Freq$GBM >0.5)] <- "High"
  print(summary(as.factor(Freq$similar)))
  C$Freq <- Freq
  high <- Freq[which(Freq$BRCA > 0.5 & Freq$GBM >0.5),]
  
  C$high <- high
  C$plot <- (qplot(BRCA,GBM, data=Freq, color=similar) +
               geom_text(aes(label=row.names(high)), position="jitter",color="black",data=high))  
  return(C)
}

#plot figures
plot121 <- Draw.overlap(1,2,1)  #mRNA
plot122 <- Draw.overlap2(1,2,2) #methylation
plot123 <- Draw.overlap(1,2,3)  #CNA

plot121$plot
plot122$plot
plot123$plot


