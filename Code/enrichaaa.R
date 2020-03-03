ENRICH<-function(i){
  x<-diff.coxgene[[i]]
  h<-as.data.frame(clusterProfiler::enrichGO(x,
                            OrgDb="org.Hs.eg.db",
                            keytype = "SYMBOL",
                            ont="BP",
                            pvalueCutoff = 1,
                            pAdjustMethod = "fdr",
                            qvalueCutoff = 1,
                            minGSSize = 20,
                            maxGSSize = 5000))
  h1<-h[h[,5]<0.005,]
  return(h1)
}

RANDOM<-function(m){
  random.ij<-list()
  set.seed(m);random.ij[[1]]<-sample(allGene,length(diff.coxgene[[i]]))
  set.seed(m+1);random.ij[[2]]<-sample(allGene,length(diff.coxgene[[j]]))
  random.enrich<-lapply(random.ij,function(x){
    h<-as.data.frame(clusterProfiler::enrichGO(x,
                              OrgDb="org.Hs.eg.db",
                              keytype = "SYMBOL",
                              ont="BP",
                              pvalueCutoff = 1,
                              pAdjustMethod = "fdr",
                              qvalueCutoff = 1,
                              minGSSize = 20,
                              maxGSSize = 5000))
    h1<-h[h[,5]<0.005,]
    return(h1)
  })
  
  result<-GOSemSim::mgoSim(random.enrich[[1]][,1],random.enrich[[2]][,1],
                    semData=d, measure="Wang")
  return(result)
}

library(doParallel)
library(foreach)
library("clusterProfiler")
library("org.Hs.eg.db")
library("GOSemSim")
load("diff_coxgene.Rdata")
load("datasets_scale.Rdata")
load("enrich_result.Rdata")
cl<- makeCluster(4)      
registerDoParallel(cl)
enrich.result<-foreach(i=1:16)%dopar%ENRICH(i)
stopCluster(cl)
save(enrich.result,file="enrich_result.Rdata")



allGene<-unique(unlist(sapply(datasets_scale,function(x){rownames(x)})))
d <- godata('org.Hs.eg.db', ont="BP", computeIC=FALSE)
tri_semsim<-c()
tri_semsim_p<-c()
cl<- makeCluster(2)      
registerDoParallel(cl)
for(i in 14)
{
  for(j  in (1:16)[-14])
  {
    sim<-mgoSim(enrich.result[[i]][,1],enrich.result[[j]][,1],semData=d, measure="Wang")
    tri_semsim<-c(tri_semsim,sim)
    random<-foreach(m=1:2)%dopar%RANDOM(m)
    random<-random[!is.na(random)]
    real.scale<-as.numeric((sim-mean(random))/sd(random))
    p<-2*pnorm(abs(real.scale),lower.tail=FALSE)
    tri_semsim_p<-c(tri_semsim_p,p)
  }
}
stopCluster(cl)
