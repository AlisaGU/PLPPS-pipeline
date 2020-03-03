bubblePlot<-function(sim,simp,title,filename){
  #sim<-semSimilirity
  #simp<-semSimilirity_p
  #title<-"title"
  #filename<-"filename.pdf"
  a<-Bubble_row_colName(sim,"all")
  a1<-data.frame(a[,1],a[,2],unlist(sim),
                 unlist(simp))
  a1[5]<-ifelse(a1[4]<0.05,"p <0.05","")
  a1<-a1[-which(is.na(a1[,3])),]
  colnames(a1)<-c("D1","D2","Similarity","p","Significance")
  Bubble(a1,title,fillStr = "fill",
         hueFirst = "#FFFFFF",hueLast = "#723C2E",
         gradientLow = "#FFD99B",gradientHigh = "#EFC000")
  ggsave(filename,width = 7,height = 6,useDingbats=F)
}