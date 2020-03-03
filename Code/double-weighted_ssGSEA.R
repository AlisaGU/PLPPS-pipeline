Project.to.GeneSet <- function(
  # data.matrix containing gene expression data
  data.array,
  
  # exponential weight applied to ranking in calculation of
  # enrichment score
  weight = 0,
  
  # weight applied to each element of functional gene set
  w.gset = 1,
  
  # gene set projecting expression data to
  gene.set)
{
  m<-data.array
  Ns<-dim(m)[2]
  Ng<-dim(m)[1]
  for (j in 1:Ns) {  # column rank normalization
    m[,j] <- rank(m[,j], ties.method = "average")
  }
  m <- 10000*m/Ng
  data.array<-m
  gene.names <- row.names(data.array)
  n.rows <- dim(data.array)[1]
  n.cols <- dim(data.array)[2]
  
  ES.vector <- vector(length=n.cols)
  ranked.expression <- vector(length=n.rows, mode="numeric")
  
  # Compute ES score for signatures in each sample
  for (sample.index in 1:n.cols) {
    # gene.list is permutation (list of row indices) of the normalized expression data, where
    # permutation places expression data in decreasing order
    # Note that in ssGSEA we rank genes by their expression level rather than by a measure of correlation
    # between expression profile and phenotype.
    gene.list <- order(data.array[, sample.index], decreasing=T)
    
    # gene.set2 contains the indices of the matching genes.
    # Note that when input GCT file is ATARiS-generated, elements of
    # gene.names may not be unique; the following code insures each element
    # of gene.names that is present in the gene.set is referenced in gene.set2
    gene.set2 <- seq(1:length(gene.names))[!is.na(match(gene.names, gene.set))]
    
    # transform the normalized expression data for a single sample into ranked (in decreasing order)
    # expression values
    if (weight == 0) {
      # don't bother doing calcuation, just set to 1
      ranked.expression <- rep(1, n.rows)
    } else if (weight > 0) {
      # calculate z.score of normalized (e.g., ranked) expression values
      x <- data.array[gene.list, sample.index]
      ranked.expression <- (x - mean(x))/sd(x)
    }
    
    # tag.indicator flags, within the ranked list of genes, those that are in the gene set
    tag.indicator <- sign(match(gene.list, gene.set2, nomatch=0))    # notice that the sign is 0 (no tag) or 1 (tag)
    no.tag.indicator <- 1 - tag.indicator
    N <- length(gene.list)
    Nh <- length(gene.set2)
    Nm <-  N - Nh
    # ind are indices into ranked.expression, whose values are in decreasing order, corresonding to
    # genes that are in the gene set
    ind = which(tag.indicator==1)
    w.gset.vec<-rep(1,N)
    # weight the ranked gene with 1 for miss or w.gset for hit
    if(length(w.gset)==length(gene.set))
    {w.gset.vec[na.omit(match(gene.set,gene.names[gene.list]))]<-w.gset[!is.na(match(gene.set,gene.names[gene.list]))]}else{print("warning!!!The default weight of gene set: 1 is used, and may the length of gene set and gene set's weight are not equal!!!")}
    ranked.expression <- (abs(ranked.expression[ind])^weight)*w.gset.vec[ind]
    
    sum.ranked.expression = sum(ranked.expression)
    # "up" represents the peaks in the mountain plot; i.e., increments in the running-sum
    up = ranked.expression/sum.ranked.expression
    # "gaps" contains the lengths of the gaps between ranked pathway genes
    gaps = (c(ind-1, N) - c(0, ind))
    # "down" contain the valleys in the mountain plot; i.e., the decrements in the running-sum
    down = gaps/Nm
    # calculate the cumulative sums at each of the ranked pathway genes
    RES = cumsum(c(up,up[Nh])-down)
    valleys = RES[1:Nh]-up
    
    max.ES = max(RES)
    min.ES = min(valleys)
    
    if( max.ES > -min.ES ){
      arg.ES <- which.max(RES)
    } else{
      arg.ES <- which.min(RES)
    }
    # calculates the area under RES by adding up areas of individual
    # rectangles + triangles
    gaps = gaps+1
    RES = c(valleys,0) * (gaps) + 0.5*( c(0,RES[1:Nh]) - c(valleys,0) ) * (gaps)
    ES = sum(RES)
    
    ES.vector[sample.index] <- ES
    
  }
  return(list(ES.vector = ES.vector))
  
} # end of Project.to.GeneSet
