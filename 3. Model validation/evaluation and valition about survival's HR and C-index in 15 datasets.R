clinicalTrain<-clinical[trainPo]
clinicalMatrix<-getSurvivalMatrix(clinicalTrain,"global","Patch et al.")
survivalMatrix<-as.data.frame(cbind(clinicalMatrix,score,1))

a<-sapply(ssgsea.out,function(x){
  np<-length(which(x[3,]>medianScore))
  nn<-length(which(x[3,]<medianScore))
  return(c(np,nn))
})
rownames(a)<-c("score>median score","score<median score")


survivalMatrix<-read.csv("survivalMatrix.csv",header = T,row.names = 1,as.is = T)
survivalMatrix_sort<-survivalMatrix[order(survivalMatrix[,3]),]
survivalMatrix_sort<-cbind(survivalMatrix_sort,
                           seq(from=1,to=nrow(survivalMatrix_sort),by=1))
colnames(survivalMatrix_sort)[c(2,5)]<-c("SurvivalStatus","order")
survivalMatrix_sort$score<-survivalMatrix_sort$score-medianScore

poL0<-which(survivalMatrix_sort$score>0)
poS0<-which(survivalMatrix_sort$score<0)
ma<-matrix(c(length(which((survivalMatrix_sort$score>0)&(survivalMatrix_sort$SurvivalStatus=="Alive"))),
             length(which((survivalMatrix_sort$score<0)&(survivalMatrix_sort$SurvivalStatus=="Alive"))),
             length(which((survivalMatrix_sort$score>0)&(survivalMatrix_sort$SurvivalStatus=="Dead"))),
             length(which((survivalMatrix_sort$score<0)&(survivalMatrix_sort$SurvivalStatus=="Dead")))),nrow = 2)
rownames(ma)<-c("High Risk","Low Risk")
colnames(ma)<-c("Alive","Dead")
chisq.test(ma)

ggplot(survivalMatrix_sort,aes(order,score))+
  geom_bar(stat = "identity",aes(colour=SurvivalStatus,fill=SurvivalStatus))+
  scale_color_manual(values=c("#F5BB2C","#033994"))+
  scale_fill_manual(values =c("#F5BB2C","#033994"))+
  theme_classic()+
  scale_y_continuous(breaks  = c(-12360,-2688.0,11330)-medianScore,
                     labels = c("0","50th","100th"))+
  labs(x="",y="Percentile of FIPI",title="Training")+
  geom_vline(xintercept = 595,linetype="dotted")+
  theme(aspect.ratio =7/10,legend.position = c(0.01,1),
        legend.justification = c(0.01,1),
        plot.title = element_text(hjust = 0.5,size=15),
        legend.text = element_text(size=12),
        legend.title = element_text(size=14),
        axis.text.x = element_text(color = "black",size=14),
        axis.text.y=element_text(color="black",size=14),
        axis.title.y=element_text(color="black",size=16))+
  annotate("text",length(score)/2-400,-6000,label="Low-Risk",color="#01939A",size=5)+
  annotate("text",length(score)/2+400,6100,label="High-Risk",color="#E73018",size=5)

ggplot(survivalMatrix_sort,aes(order,time))+
  geom_point(aes(color=SurvivalStatus))+
  labs(x="",y="Survival Time (Day)")+
  geom_vline(xintercept = 1686/2)+
  theme_classic()+
  theme(aspect.ratio = 9/10,
        legend.position = c(0.01,1),
        legend.justification = c(0.01,1),
        legend.text = element_text(size=12),
        legend.title = element_text(size=14),
        axis.text.x = element_text(color = "black",size=14),
        axis.text.y=element_text(color="black",size=14),
        axis.title.y=element_text(color="black",size=16))+
  scale_colour_manual(values =c("#FFAC40",rgb(0,0,0)))

Data<-getHeatmapExp(gene,datasets_nonscale[trainPo])
Data1<-Data[,rownames(survivalMatrix_sort)]
annotation_row<-data.frame(Gene=c(rep("Protect",25),rep("Risk",43)))
rownames(annotation_row)<-gene
annotation_col<-data.frame(Sample=c(rep("High Risk",595),
                                    rep("Low Risk",595)))
rownames(annotation_col)<-colnames(Data1)
pheatmap(Data1,cluster_rows = F,cluster_cols = F,scale="row",
         show_colnames = F,annotation_row = annotation_row,
         annotation_col = annotation_col,na_col="green")
         
clinicalTest<-clinical[testPo]
clinicalMatrix<-getSurvivalMatrix(clinicalTest,"global","Patch et al.")
survivalMatrix<-data.frame(time=as.numeric(clinicalMatrix[,1]),
                           status=as.character(clinicalMatrix[,2]),
                           score=as.numeric(score.test))
survivalMatrix_sort<-survivalMatrix[order(survivalMatrix[,3]),]
survivalMatrix_sort<-cbind(survivalMatrix_sort,
                           seq(from=1,to=nrow(survivalMatrix_sort),by=1))
colnames(survivalMatrix_sort)[c(2,4)]<-c("SurvivalStatus","order")
survivalMatrix_sort$score<-survivalMatrix_sort$score-medianScore#medianScore

poL0<-which(survivalMatrix_sort$score>0)
poS0<-which(survivalMatrix_sort$score<0)
ma<-matrix(c(length(which((survivalMatrix_sort$score>0)&(survivalMatrix_sort$SurvivalStatus=="Alive"))),
             length(which((survivalMatrix_sort$score<0)&(survivalMatrix_sort$SurvivalStatus=="Alive"))),
             length(which((survivalMatrix_sort$score>0)&(survivalMatrix_sort$SurvivalStatus=="Dead"))),
             length(which((survivalMatrix_sort$score<0)&(survivalMatrix_sort$SurvivalStatus=="Dead")))),nrow = 2)
rownames(ma)<-c("High Risk","Low Risk")
colnames(ma)<-c("Alive","Dead")
chisq.test(ma)


ggplot(survivalMatrix_sort,aes(order,score))+
  geom_bar(stat = "identity",aes(colour=SurvivalStatus,fill=SurvivalStatus))+
  scale_color_manual(values=c("#F5BB2C","#033994"))+
  scale_fill_manual(values =c("#F5BB2C","#033994"))+
  theme_classic()+
  labs(x="",y="",title="Test")+
  geom_vline(xintercept = 195,linetype="dotted")+
  theme(aspect.ratio =7/10,legend.position = c(0.01,1),
        legend.justification = c(0.01,1),
        plot.title = element_text(hjust = 0.5,size=15),
        legend.text = element_text(size=12),
        legend.title = element_text(size=14),
        axis.text.x = element_text(color = "black",size=14),
        #axis.text.y=element_blank(),
        axis.title.y=element_text(color="black",size=16))+
  annotate("text",100,-6000,label="Low-Risk",color="#01939A",size=5)+
  annotate("text",350,6100,label="High-Risk",color="#E73018",size=5)
  
clinicalIndepend<-clinical[independPo]#TCGA、patch
clinicalMatrix11<-getSurvivalMatrix1(clinicalIndepend[1],"global","Patch et al.")
clinicalMatrix22<-getSurvivalMatrix(clinicalIndepend[2],"global","Patch et al.")
clinicalMatrix<-rbind(clinicalMatrix11,clinicalMatrix22)
group.test.independ<-c()
for(i in 1:length(score.test.independ)){
  group.test.i<-group(score.test.independ[i],type="0",medianScore)
  group.test.independ<-c(group.test.independ,group.test.i)
}
survivalMatrix<-data.frame(time=as.numeric(clinicalMatrix[,1]),
                           status=as.character(clinicalMatrix[,2]),
                           score=as.numeric(score.test.independ),
                           group=group.test.independ)
names(survivalMatrix)<-c("time","status","score","group")
survivalMatrix$group<-ifelse(survivalMatrix$group=="negative","High Risk","Low Risk")
survivalMatrix_sort<-survivalMatrix[order(survivalMatrix[,3]),]
survivalMatrix_sort<-cbind(survivalMatrix_sort,
                           seq(from=1,to=nrow(survivalMatrix_sort),by=1))
colnames(survivalMatrix_sort)[c(2,5)]<-c("SurvivalStatus","order")
survivalMatrix_sort$SurvivalStatus<-factor(survivalMatrix_sort$SurvivalStatus)
survivalMatrix_sort$score<-survivalMatrix_sort$score-medianScore
ggplot(survivalMatrix_sort,aes(order,score))+
  geom_bar(stat = "identity",aes(colour=SurvivalStatus,fill=SurvivalStatus))+
  scale_color_manual(values=c("#F5BB2C","#033994"))+
  scale_fill_manual(values =c("#F5BB2C","#033994"))+
  theme_classic()+
  labs(x="",y="",title="Independent Test")+
  geom_vline(xintercept = 159,linetype="dotted")+
  theme(aspect.ratio =7/10,legend.position = c(0.01,1),
        legend.justification = c(0.01,1),
        legend.text = element_text(size=12),
        plot.title = element_text(hjust = 0.5,size=15),
        legend.title = element_text(size=14),
        axis.text.x = element_text(color = "black",size=14),
        axis.text.y=element_blank(),
        axis.title.y=element_text(color="black",size=16))+
  annotate("text",100,-4000,label="Low-Risk",color="#01939A",size=5)+
  annotate("text",250,4000,label="High-Risk",color="#E73018",size=5)
  
  
  
for(i in 1:nrow(survivalMatrix)){
  survivalMatrix[i,4]<-group(survivalMatrix[i,3],type="0",medianScore)
}
names(survivalMatrix)<-c("time","status","score","group")
survivalMatrix$status<-ifelse(survivalMatrix$status=="Alive",0,1)
survivalMatrix$group<-ifelse(survivalMatrix$group=="negative","High Risk","Low Risk")
eee<-table(survivalMatrix$status,survivalMatrix$group)
chisq.test(eee)

fit<-survfit(Surv(time,status)~group,data=survivalMatrix)
diff<-survdiff(Surv(time,status)~group,data=survivalMatrix)

ggsurvplot(fit, data = survivalMatrix, title="Training dataset(n = 1190)",
           risk.table = TRUE,risk.table.col = "strata",pval = 1e-10,
           risk.table.height = 0.2,xlab="Time(Days)",legend.title="",
           legend.labs=c("High Risk","Low Risk"),risk.table.title="No. at risk",
           palette = c("#E7201C","#0A499E"),
           ggtheme = theme_bw()+theme(panel.border = element_blank(),panel.grid=element_blank(),
                                      axis.line = element_line(colour = "black"),
                                      axis.text.x = element_text(size=15,color = "black"),
                                      axis.text.y=element_text(size=15,color="black"),
                                      plot.title = element_text(colour = "black", face = "bold", 
                                                                size = 14, vjust = 1,hjust = 0.5),
                                      axis.title.x = element_text(size=15,color = "black"),
                                      axis.title.y = element_text(size=15,color = "black"),
                                      legend.text = element_text(size = 15),
                                      legend.title = element_text(size = 15)))
                                      
                                      
survivalMatrix<-read.csv("survivalMatrix.csv",header = T,row.names = 1,as.is = T)
for(i in 1:nrow(survivalMatrix)){
  survivalMatrix[i,4]<-groupQ(survivalMatrix[i,3],Q1score,medianScore,Q3score)
}
names(survivalMatrix)<-c("time","status","score","group")
survivalMatrix$status<-ifelse(survivalMatrix$status=="Alive",0,1)
survivalMatrix$group<-factor(survivalMatrix$group,levels = c("Class one","Class two",
                                                             "Class three","Class four"))
eee<-table(survivalMatrix$status,survivalMatrix$group)
chisq.test(eee)
fit<-survfit(Surv(time,status)~group,data=survivalMatrix)
diff<-survdiff(Surv(time,status)~group,data=survivalMatrix)#2.22e-16

ggsurvplot(fit, data = survivalMatrix, title="Training dataset(n = 1190)",
           risk.table = TRUE,risk.table.col = "strata",pval = T,
           risk.table.height = 0.2,xlab="Time(Days)",legend.title="",
           risk.table.title="No. at risk",
           ggtheme = theme_bw()+theme(panel.border = element_blank(),panel.grid=element_blank(),
                                      axis.line = element_line(colour = "black"),
                                      axis.text.x = element_text(size=15,color = "black"),
                                      axis.text.y=element_text(size=15,color="black"),
                                      plot.title = element_text(colour = "black", face = "bold", 
                                                                size = 14, vjust = 1,hjust = 0.5),
                                      axis.title.x = element_text(size=15,color = "black"),
                                      axis.title.y = element_text(size=15,color = "black"),
                                      legend.text = element_text(size = 15),
                                      legend.title = element_text(size = 15)))
                                      
                                      
clinicalTrain<-clinical[trainPo]
for(i in 1:length(ssgsea.out)){
  survivalMatrix.i<-getSurvivalMatrix(clinicalTrain[[i]],
                                      "seperate","Patch et al.",
                                      names(clinicalTrain)[i])
  survivalMatrix.ii<-cbind(survivalMatrix.i,ssgsea.out[[i]][3,],1)
  for(j in 1:nrow(survivalMatrix.ii)){
    survivalMatrix.ii[j,4]<-group(as.numeric(survivalMatrix.ii[j,3]),
                                  type="0",medianScore)
  }
  survivalMatrix.iii<-data.frame(time=as.numeric(survivalMatrix.ii[,1]),
                                 status=as.character(survivalMatrix.ii[,2]),
                                 score=as.numeric(survivalMatrix.ii[,3]),
                                 group=as.character(survivalMatrix.ii[,4]))
  survivalMatrix.iii$status<-ifelse(survivalMatrix.iii$status=="Alive",0,1)
  survivalMatrix.iii$group<-ifelse(survivalMatrix.iii$group=="negative",
                                   "Low Risk","High Risk")
  fit<-survfit(Surv(time,as.numeric(status))~group,
               data=survivalMatrix.iii)
  diff<-survdiff(Surv(time,as.numeric(status))~group,
                 data=survivalMatrix.iii)
  print(names(ssgsea.out)[i])
  print(fit)
  print(diff)
  ggsurvplot(fit, data = survivalMatrix.iii, title=paste(names(ssgsea.out)[i],"(n = ",
                                                         sapply(ssgsea.out,ncol)[i],")",sep=""),
             risk.table = TRUE,risk.table.col = "strata",pval = T,
             risk.table.height = 0.2,xlab="Time(Days)",legend.title="",
             legend.labs=c("High Risk","Low Risk"),risk.table.title="No. at risk",
             palette = c("#E7201C","#0A499E"),
             ggtheme = theme_bw()+theme(panel.border = element_blank(),panel.grid=element_blank(),
                                        axis.line = element_line(colour = "black"),
                                        axis.text.x = element_text(size=15,color = "black"),
                                        axis.text.y=element_text(size=15,color="black"),
                                        plot.title = element_text(colour = "black", face = "bold", 
                                                                  size = 14, vjust = 1,hjust = 0.5),
                                        axis.title.x = element_text(size=15,color = "black"),
                                        axis.title.y = element_text(size=15,color = "black"),
                                        legend.text = element_text(size = 15),
                                        legend.title = element_text(size = 15)))
}

clinicalMatrix<-getSurvivalMatrix(clinical[14],"global","Patch et al.")
colnames(clinicalMatrix)<-c("time","status")
for(i in 1:nrow(clinicalMatrix)){
  if(clinicalMatrix[i,2]=="Dead")
    clinicalMatrix[i,2]<-1
  if(clinicalMatrix[i,2]=="Alive")
    clinicalMatrix[i,2]<-"0"
}
group.test<-c()
for(i in 1:length(score.test.tcga.cel)){
  group.test.i<-group(score.test.tcga.cel[i],type="0",medianScore)
  group.test<-c(group.test,group.test.i)
}
survivalMatrix<-data.frame(time=as.numeric(clinicalMatrix[,1]),
                           status=as.character(clinicalMatrix[,2]),
                           score=score.test.tcga.cel,
                           group=group.test)
names(survivalMatrix)<-c("time","status","score","group")
survivalMatrix$group<-ifelse(survivalMatrix$group=="negative","High Risk","Low Risk")
summary(coxph(Surv(time,as.numeric(status))~ score,survivalMatrix))$coefficients[-4]
fit<-survfit(Surv(time,as.numeric(status))~group,data=survivalMatrix)
diff<-survdiff(Surv(time,as.numeric(status))~group,data=survivalMatrix)
ggsurvplot(fit, data = survivalMatrix, title="Test dataset(n = 495)",
           risk.table = TRUE,risk.table.col = "strata",
           pval =  0.0343,
           risk.table.height = 0.2,xlab="Time(Days)",legend.title="",
           legend.labs=c("High Risk","Low Risk"),risk.table.title="No. at risk",
           palette = c("#E7201C","#0A499E"),
           ggtheme = theme_bw()+theme(panel.border = element_blank(),panel.grid=element_blank(),
                                      axis.line = element_line(colour = "black"),
                                      axis.text.x = element_text(size=15,color = "black"),
                                      axis.text.y=element_text(size=15,color="black"),
                                      plot.title = element_text(colour = "black", face = "bold", 
                                                                size = 14, vjust = 1,hjust = 0.5),
                                      axis.title.x = element_text(size=15,color = "black"),
                                      axis.title.y = element_text(size=15,color = "black"),
                                      legend.text = element_text(size = 15),
                                      legend.title = element_text(size = 15)))
                                      
                                      
clinicalTest<-clinical[testPo]
clinicalMatrix<-getSurvivalMatrix(clinicalTest,"global","Patch et al.")
group.test<-c()
for(i in 1:length(score.test)){
  group.test.i<-group(score.test[i],type="0",medianScore)
  group.test<-c(group.test,group.test.i)
}
survivalMatrix<-data.frame(time=as.numeric(clinicalMatrix[,1]),
                           status=as.character(clinicalMatrix[,2]),
                           score=as.numeric(score.test),
                           group=group.test)
names(survivalMatrix)<-c("time","status","score","group")
survivalMatrix$status<-ifelse(survivalMatrix$status=="Alive",0,1)
survivalMatrix$group<-ifelse(survivalMatrix$group=="negative","High Risk","Low Risk")
eee<-table(survivalMatrix$status,survivalMatrix$group)
chisq.test(eee)
fit<-survfit(Surv(time,as.numeric(status))~group,data=survivalMatrix)
diff<-survdiff(Surv(time,as.numeric(status))~group,data=survivalMatrix)

ggsurvplot(fit, data = survivalMatrix, title="Test dataset(n = 438)",
           risk.table = TRUE,risk.table.col = "strata",
           pval =  0.000364,
           risk.table.height = 0.2,xlab="Time(Days)",legend.title="",
           legend.labs=c("High Risk","Low Risk"),risk.table.title="No. at risk",
           palette = c("#E7201C","#0A499E"),
           ggtheme = theme_bw()+theme(panel.border = element_blank(),panel.grid=element_blank(),
                                      axis.line = element_line(colour = "black"),
                                      axis.text.x = element_text(size=15,color = "black"),
                                      axis.text.y=element_text(size=15,color="black"),
                                      plot.title = element_text(colour = "black", face = "bold", 
                                                                size = 14, vjust = 1,hjust = 0.5),
                                      axis.title.x = element_text(size=15,color = "black"),
                                      axis.title.y = element_text(size=15,color = "black"),
                                      legend.text = element_text(size = 15),
                                      legend.title = element_text(size = 15)))
                                      
clinicalTest<-clinical[testPo]
clinicalMatrix<-getSurvivalMatrix(clinicalTest,"global","Patch et al.")
group.test<-c()
for(i in 1:length(score.test)){
  group.test.i<-groupQ(score.test[i],Q1score,medianScore,Q3score)
  group.test<-c(group.test,group.test.i)
}
survivalMatrix<-data.frame(time=as.numeric(clinicalMatrix[,1]),
                           status=as.character(clinicalMatrix[,2]),
                           score=as.numeric(score.test),
                           group=group.test)
names(survivalMatrix)<-c("time","status","score","group")
survivalMatrix$status<-ifelse(survivalMatrix$status=="Alive",0,1)
survivalMatrix$group<-factor(survivalMatrix$group,levels = c("Class one","Class two",
                                                             "Class three","Class four"))
eee<-table(survivalMatrix$status,survivalMatrix$group)
chisq.test(eee)
fit<-survfit(Surv(time,as.numeric(status))~group,data=survivalMatrix)
diff<-survdiff(Surv(time,as.numeric(status))~group,data=survivalMatrix)

ggsurvplot(fit, data = survivalMatrix, title="Test dataset(n = 438)",
           risk.table = TRUE,risk.table.col = "strata",
           pval =  T,
           risk.table.height = 0.2,xlab="Time(Days)",legend.title="",
           risk.table.title="No. at risk",
           ggtheme = theme_bw()+theme(panel.border = element_blank(),panel.grid=element_blank(),
                                      axis.line = element_line(colour = "black"),
                                      axis.text.x = element_text(size=15,color = "black"),
                                      axis.text.y=element_text(size=15,color="black"),
                                      plot.title = element_text(colour = "black", face = "bold", 
                                                                size = 14, vjust = 1,hjust = 0.5),
                                      axis.title.x = element_text(size=15,color = "black"),
                                      axis.title.y = element_text(size=15,color = "black"),
                                      legend.text = element_text(size = 15),
                                      legend.title = element_text(size = 15)))
                                      
                                      
clinicalIndepend<-clinical[independPo]#TCGA、patch
clinicalMatrix11<-getSurvivalMatrix1(clinicalIndepend[1],"global","Patch et al.")
clinicalMatrix22<-getSurvivalMatrix(clinicalIndepend[2],"global","Patch et al.")
clinicalMatrix<-rbind(clinicalMatrix11,clinicalMatrix22)
group.test.independ<-c()
for(i in 1:length(score.test.independ)){
  group.test.i<-group(score.test.independ[i],type="0",medianScore)
  group.test.independ<-c(group.test.independ,group.test.i)
}
survivalMatrix<-data.frame(time=as.numeric(clinicalMatrix[,1]),
                           status=as.character(clinicalMatrix[,2]),
                           score=as.numeric(score.test.independ),
                           group=group.test.independ)
names(survivalMatrix)<-c("time","status","score","group")
survivalMatrix$group<-ifelse(survivalMatrix$group=="negative","High Risk","Low Risk")
eee<-table(survivalMatrix$status,survivalMatrix$group)
chisq.test(eee)
fit<-survfit(Surv(time,as.numeric(status))~group,data=survivalMatrix)
diff<-survdiff(Surv(time,as.numeric(status))~group,data=survivalMatrix)

ggsurvplot(fit, data = survivalMatrix, title=" Independent Test Dataset(n = 354)",
           risk.table = TRUE,risk.table.col = "strata",
           pval =  0.00359,
           risk.table.height = 0.2,xlab="Time(Days)",legend.title="",
           legend.labs=c("High Risk","Low Risk"),risk.table.title="No. at risk",
           palette = c("#E7201C","#0A499E"),
           ggtheme = theme_bw()+theme(panel.border = element_blank(),panel.grid=element_blank(),
                                      axis.line = element_line(colour = "black"),
                                      axis.text.x = element_text(size=15,color = "black"),
                                      axis.text.y=element_text(size=15,color="black"),
                                      plot.title = element_text(colour = "black", face = "bold", 
                                                                size = 14, vjust = 1,hjust = 0.5),
                                      axis.title.x = element_text(size=15,color = "black"),
                                      axis.title.y = element_text(size=15,color = "black"),
                                      legend.text = element_text(size = 15),
                                      legend.title = element_text(size = 15)))
dev.off()

#四分
clinicalIndepend<-clinical[independPo]#TCGA、patch
clinicalMatrix11<-getSurvivalMatrix1(clinicalIndepend[1],"global","Patch et al.")
clinicalMatrix22<-getSurvivalMatrix(clinicalIndepend[2],"global","Patch et al.")
clinicalMatrix<-rbind(clinicalMatrix11,clinicalMatrix22)
group.test.independ<-c()
for(i in 1:length(score.test.independ)){
  group.test.i<-groupQ(score.test.independ[i],Q1score,medianScore,Q3score)
  group.test.independ<-c(group.test.independ,group.test.i)
}
survivalMatrix<-data.frame(time=as.numeric(clinicalMatrix[,1]),
                           status=as.character(clinicalMatrix[,2]),
                           score=as.numeric(score.test.independ),
                           group=group.test.independ)
names(survivalMatrix)<-c("time","status","score","group")
survivalMatrix$group<-factor(survivalMatrix$group,levels = c("Class one","Class two",
                                                             "Class three","Class four"))
eee<-table(survivalMatrix$status,survivalMatrix$group)
chisq.test(eee)
fit<-survfit(Surv(time,as.numeric(status))~group,data=survivalMatrix)
diff<-survdiff(Surv(time,as.numeric(status))~group,data=survivalMatrix)

ggsurvplot(fit, data = survivalMatrix, title=" Independent Test Dataset(n = 354)",
           risk.table = TRUE,risk.table.col = "strata",
           pval = T,
           risk.table.height = 0.2,xlab="Time(Days)",legend.title="",
           risk.table.title="No. at risk",
           ggtheme = theme_bw()+theme(panel.border = element_blank(),panel.grid=element_blank(),
                                      axis.line = element_line(colour = "black"),
                                      axis.text.x = element_text(size=15,color = "black"),
                                      axis.text.y=element_text(size=15,color="black"),
                                      plot.title = element_text(colour = "black", face = "bold", 
                                                                size = 14, vjust = 1,hjust = 0.5),
                                      axis.title.x = element_text(size=15,color = "black"),
                                      axis.title.y = element_text(size=15,color = "black"),
                                      legend.text = element_text(size = 15),
                                      legend.title = element_text(size = 15)))
                                      
clinicalTest<-clinical[-trainPo]

for(i in 1:length(ssgsea.out.test)){
  survivalMatrix.i<-getSurvivalMatrix(clinicalTest[[i]],
                                      "seperate","Patch et al.",
                                      names(clinicalTest)[i])
  survivalMatrix.ii<-cbind(survivalMatrix.i,ssgsea.out.test[[i]][4,],1)
  for(j in 1:nrow(survivalMatrix.ii)){
    survivalMatrix.ii[j,4]<-group(as.numeric(survivalMatrix.ii[j,3]),
                                  type="0",summary(score)[-4])
  }
  survivalMatrix.iii<-data.frame(time=as.numeric(survivalMatrix.ii[,1]),
                                 status=as.character(survivalMatrix.ii[,2]),
                                 score=as.numeric(survivalMatrix.ii[,3]),
                                 group=as.character(survivalMatrix.ii[,4]))
  fit<-survfit(Surv(time,as.numeric(status))~group,
               data=survivalMatrix.iii)
  ggsurvplot(fit, data = survivalMatrix.iii, size = 1,
             risk.table = TRUE,risk.table.col = "strata",pval = TRUE,
             risk.table.height = 0.2,
             ggtheme = theme_bw())
}

sccForestData<-rbind(sstselfEval("trainPo","Training"),
                     sstselfEval("testPo","Test"),sstselfEval("independent","Independent"))
rownames(sccForestData)<-paste(c(rep("Training",5),rep("Test",5),rep("Independent",5)),
                               rownames(sccForestData),sep="-")
sccForestData<-as.data.frame(sccForestData)
Data<-data.frame(lnHR=sccForestData$lnHR,
                 selnHR=sccForestData$selnHR)
rownames(Data)<-rownames(sccForestData)
pdf("sccForest.pdf",width=6,height=4)
forest(metagen(Data$lnHR,Data$selnHR,sm="HR"),overall = F,squaresize = 0.4,
       weight.study = "same",xlim = c(0.5,4),leftcols = "studlab",col.square="#2D5662",
       studlab=rownames(Data),rightcols = c("effect.ci"),ref=1,plotwidth = "3cm",
       col.predict.lines="#2D5662")
dev.off()
HR<-sccForestData$HR
lnHR<-sccForestData$lnHR
selnHR<-sccForestData$selnHR
low.HR<-exp(lnHR-1.96*selnHR)
high.HR<-exp(lnHR+1.96*selnHR)
p<-sccForestData$`p value`
a<-data.frame(sccForestData$`Low Risk`,
              sccForestData$`High Risk`,
              format(HR,digits=2),
              paste(format(low.HR,digits=2),
                    format(high.HR,digits=2),sep="-"),
              p=format(p,digits=2,scientific = 2))
colnames(a)<-c("Low Risk","High Risk","HR","95%CI of HR","P value")
rownames(a)<-rownames(sccForestData)



ooo<-sstselfEval("global","global")
sccForestData<-ooo[-1,]
sccForestData<-as.data.frame(sccForestData)
HR<-sccForestData$HR
lnHR<-sccForestData$lnHR
selnHR<-sccForestData$selnHR
low.HR<-exp(lnHR-1.96*selnHR)
high.HR<-exp(lnHR+1.96*selnHR)
p<-sccForestData$`p value`
a<-data.frame(sccForestData$`Low Risk`,
              sccForestData$`High Risk`,
              round(HR,digits=2),
              paste(round(low.HR,digits=2),
                    round(high.HR,digits=2),sep="-"),
              p=format(p,digits=2,scientific = 2))
colnames(a)<-c("Low Risk","High Risk","HR","95%CI of HR","P value")
rownames(a)<-rownames(sccForestData)



sccForestData<-rbind(sstselfEval("trainPo","Training"),
                     sstselfEval("TestALL","TestALL"))
                     
unicox<-rbind(sstselfEvalUM("trainPo")[[1]],
              sstselfEvalUM("testPo")[[1]],
              sstselfEvalUM("independent")[[1]])
rownames(unicox)<-paste(c(rep("Training",3),rep("Test",3),rep("Independent",3)),
                        rownames(unicox),sep="-")
                        
multicox<-rbind(sstselfEvalUM("trainPo")[[2]],
                sstselfEvalUM("testPo")[[2]],
                sstselfEvalUM("independent")[[2]])
rownames(multicox)<-paste(c(rep("Training",3),rep("Test",3),rep("Independent",3)),
                          rownames(multicox),sep="-")
                          
multicox<-as.data.frame(multicox)
Data<-data.frame(lnHR=multicox$coef,
                 selnHR=multicox$`se(coef)`)
rownames(Data)<-rownames(multicox)
forest(metagen(Data$lnHR,Data$selnHR,sm="HR"),overall = F,squaresize = 0.4,
       weight.study = "same",xlim = c(0.5,2.5),leftcols = "studlab",col.square="#2D5662",
       studlab=rownames(Data),rightcols = c("effect.ci"),ref=1,plotwidth = "3cm",
       col.predict.lines="#2D5662")
       
       
unicox<-rbind(sstselfEvalUMContinue("trainPo")[[1]],
              sstselfEvalUMContinue("testPo")[[1]],
              sstselfEvalUMContinue("independent")[[1]])
rownames(unicox)<-paste(c(rep("Training",3),rep("Test",3),rep("Independent",3)),
                        rownames(unicox),sep="-")
                        
gradeStatus<-c(0,rep(1,5),0,rep(1,8),2)
clinicalTrain<-clinical[trainPo]
gradeStatusTrain<-gradeStatus[trainPo]
result<-HGClinicalStatus(clinicalTrain,gradeStatusTrain)
idx.tumorstage.train<-result[,1]
idx.tumorgrade.train<-result[,2]
idx.age.train<-result[,3]
stageNameTrain<-c(rep("tumorstage",15),"tumour_stage")[trainPo]
gradeNameTrain<-c(rep("grade",15),"tumour_grade")[trainPo]
ageNameTrain<-c(rep("age_at_initial_pathologic_diagnosis",15),
                "donor_age_at_diagnosis")[trainPo]
timeNameTrain<-c(rep("days_to_death",15),"donor_survival_time")[trainPo]
eventNameTrain<-c(rep("vital_status",15),"donor_vital_status")[trainPo]
selfEvalResult<-list()
par(mfrow=c(4,4))
for(i in 1:length(clinicalTrain)){#length(clinicalTrain)
  score.i<-scale(as.numeric(ssgsea.out[[i]][4,]))
  phe<-clinicalTrain[[i]]
  
  sCol<-which(colnames(phe)==stageNameTrain[i])
  gCol<-which(colnames(phe)==gradeNameTrain[i])
  aCol<-which(colnames(phe)==ageNameTrain[i])
  timeCol<-which(colnames(phe)==timeNameTrain[i])
  eventCol<-which(colnames(phe)==eventNameTrain[i])
  
  result<-selfEval(score.i,phe,names(clinicalTrain)[i],
                   idx.tumorstage.train[i],idx.tumorgrade.train[i],
                   idx.age.train[i],sCol,gCol,aCol,timeCol,eventCol)
  selfEvalResult[[names(clinicalTrain)[i]]]<-result
}
selfEvalCoxResult<-c()
for(i in 1:length(selfEvalResult)){
  x<-selfEvalResult[[i]]
  Xname<-names(selfEvalResult)[i]
  selfEvalCoxResult<-rbind(selfEvalCoxResult,combineCoxResult(x,Xname))
}
threeYearSurvival<-sapply(selfEvalResult,function(x){x[[3]][1]})
fiveYearSurvival<-sapply(selfEvalResult,function(x){x[[3]][2]})
c.index<-sapply(selfEvalResult,function(x){x[[3]][3]})
sd.c.index<-sapply(selfEvalResult,function(x){x[[3]][4]})
low.ci<-c.index-1.96*sd.c.index
high.ci<-c.index+1.96*sd.c.index

roc.ci<-data.frame(names(selfEvalResult),format(threeYearSurvival,digits = 4),
                   format(fiveYearSurvival,digits = 4),
                   format(c.index,digits = 4),
                   paste(format(low.ci,digits = 4),
                         format(high.ci,digits = 4),sep="-"))
names(roc.ci)<-c("Datasets","Three-year survival",
                 "Five-year survival","C-index","95%CI of C-index")

for(i in 1:length(selfEvalResult)){
  a<-selfEvalResult[[i]][[1]]
  rownames(a)[grepl("sCo",rownames(a),fixed = T)]<-"stage"
  rownames(a)[grepl("aCol",rownames(a),fixed = T)]<-"age"
  rownames(a)[grepl("gCol",rownames(a),fixed = T)]<-"grade"
  a[,5]<-round(a[,5],digits = 4)
  colnames(a)<-c("lnhr","hr","se_lnhr","z","p")
  a<-as.data.frame(a)
  meta.i<-metagen(lnhr,se_lnhr,sm="HR",data = a,studlab = rownames(a))
  forest(meta.i, leftcols = c("studlab"),leftlabs = c("Study Element"),
         rightcols = c("effect.ci","p"),fontsize = 13,
         rightlabs = "p",just.studlab = "center",
         overall = F,col.square = "#363886",lwd=2,
         weight.study = "same",col.inside="black",xlab=names(selfEvalResult)[i])
}
HR<-sapply(selfEvalResult,function(x){x[[2]][1,2]})
lnHR<-sapply(selfEvalResult,function(x){x[[2]][1,1]})
selnHR<-sapply(selfEvalResult,function(x){x[[2]][1,3]})
low.HR<-exp(lnHR-1.96*selnHR)
high.HR<-exp(lnHR+1.96*selnHR)
C_index<-sapply(selfEvalResult,function(x){x[[3]][3]})
se_C_index<-sapply(selfEvalResult,function(x){x[[3]][4]})
low_c_index<-C_index-1.96*se_C_index
high_c_index<-C_index+1.96*se_C_index
HR_C<-data.frame(HR,lnHR,selnHR,low.HR,high.HR,
                 C_index,se_C_index,low_c_index,high_c_index)

lnHR<-meta.summaries(HR_C[,2],HR_C[,3],method="fixed",logscale=T)
if(lnHR$het[3]>0.05){
  summary_HR<-c(lnHR$summary,lnHR$se.summary)
}else{
  lnHR<-meta.summaries(HR_C[,2],HR_C[,3],method="random",logscale=T)
  summary_HR<-c(lnHR$summary,lnHR$se.summary)
}
C<-meta.summaries(HR_C[,6],HR_C[,7],method="fixed",logscale=F)
if(C$het[3]>0.05){
  summary_C<-c(C$summary,C$se.summary)
}else{
  C<-meta.summaries(HR_C[,6],HR_C[,7],method="random",logscale=F)
  summary_C<-c(C$summary,C$se.summary)
}
ssgsea_colla_Train<-c(summary_HR,summary_C)

ol<-colorRampPalette(c("#363886","#9ECAD4"))
a<-metagen(HR_C$lnHR,HR_C$selnHR,sm="HR")
data<-data.frame(HR=c(HR_C$lnHR,a$TE.random),
                 seHR=c(HR_C$selnHR,a$seTE.random))
rownames(data)<-c(rownames(HR_C),"Total")
forest(metagen(data$HR,data$seHR,sm="HR"),overall = F,squaresize = 0.4,
       weight.study = "same",xlim = c(1,2.5),leftcols = "studlab",col.square="#363886",
       studlab=rownames(data),rightcols = c("effect.ci"),ref=1,plotwidth = "3cm",
       at=c(1,1.5,2,2.5))
       
a<-metagen(HR_C$lnHR,HR_C$selnHR,sm="lnHR")
data<-data.frame(lnHR=c(HR_C$lnHR,a$TE.random),
                 selnHR=c(HR_C$selnHR,a$seTE.random))
rownames(data)<-c(rownames(HR_C),"Total")
forest(metagen(data$lnHR,data$selnHR,sm="lnHR"),overall = F,squaresize = 0.4,
       weight.study = "same",xlim = c(0,1),leftcols = "studlab",col.square="#363886",
       studlab=rownames(data),rightcols = c("effect.ci"),at=c(0,0.2,0.4,0.6,0.8,1),plotwidth="3.5cm")


a<-metagen(HR_C$C_index,HR_C$se_C_index,sm="C-index")
data<-data.frame(C_index=c(HR_C$C_index,a$TE.random),
                 se_C_index=c(HR_C$se_C_index,a$seTE.random))
rownames(data)<-c(rownames(HR_C),"Total")
forest(metagen(data$C_index,data$se_C_index,sm="C-index"),overall = F,squaresize = 0.4,
       weight.study = "same",xlim = c(0.5,0.8),leftcols = "studlab",col.square="#363886",
       studlab=rownames(data),rightcols = c("effect.ci"),ref=0.5,plotwidth = "3cm",
       at=c(0.5,0.6,0.7,0.8))
       
       
gradeStatus<-c(0,rep(1,5),0,rep(1,8),2)
clinicalTest<-clinical[-trainPo]
gradeStatusTest<-gradeStatus[-trainPo]
result<-HGClinicalStatus(clinicalTest,gradeStatusTest)
idx.tumorstage.test<-result[,1]
idx.tumorgrade.test<-result[,2]
idx.age.test<-result[,3]
stageNameTest<-c(rep("tumorstage",15),"tumour_stage")[-trainPo]
gradeNameTest<-c(rep("grade",15),"tumour_grade")[-trainPo]
ageNameTest<-c(rep("age_at_initial_pathologic_diagnosis",15),
               "donor_age_at_diagnosis")[-trainPo]
timeNameTest<-c(rep("days_to_death",15),"donor_survival_time")[-trainPo]
eventNameTest<-c(rep("vital_status",15),"donor_vital_status")[-trainPo]
selfEvalResult<-list()
par(mfrow=c(2,2))#要手动存的哟,名字：model ROC Test
for(i in 1:length(clinicalTest)){
  score.i<-scale(as.numeric(ssgsea.out.test[[i]][4,]))
  phe<-clinicalTest[[i]]
  
  sCol<-which(colnames(phe)==stageNameTest[i])
  gCol<-which(colnames(phe)==gradeNameTest[i])
  aCol<-which(colnames(phe)==ageNameTest[i])
  timeCol<-which(colnames(phe)==timeNameTest[i])
  eventCol<-which(colnames(phe)==eventNameTest[i])
  
  result<-selfEval(score.i,phe,names(clinicalTest)[i],
                   idx.tumorstage.test[i],idx.tumorgrade.test[i],
                   idx.age.test[i],sCol,gCol,aCol,timeCol,eventCol)
  selfEvalResult[[names(clinicalTest)[i]]]<-result
}

selfEvalCoxResult<-c()
for(i in 1:length(selfEvalResult)){
  x<-selfEvalResult[[i]]
  Xname<-names(selfEvalResult)[i]
  selfEvalCoxResult<-rbind(selfEvalCoxResult,combineCoxResult(x,Xname))
}
threeYearSurvival<-sapply(selfEvalResult,function(x){x[[3]][1]})
fiveYearSurvival<-sapply(selfEvalResult,function(x){x[[3]][2]})
c.index<-sapply(selfEvalResult,function(x){x[[3]][3]})
sd.c.index<-sapply(selfEvalResult,function(x){x[[3]][4]})
low.ci<-c.index-1.96*sd.c.index
high.ci<-c.index+1.96*sd.c.index

roc.ci<-data.frame(names(selfEvalResult),format(threeYearSurvival,digits = 4),
                   format(fiveYearSurvival,digits = 4),
                   format(c.index,digits = 4),
                   paste(format(low.ci,digits = 4),
                         format(high.ci,digits = 4),sep="-"))
names(roc.ci)<-c("Datasets","Three-year survival",
                 "Five-year survival","C-index","95%CI of C-index")
                 
                 
for(i in 1:length(selfEvalResult)){
  a<-selfEvalResult[[i]][[1]]
  rownames(a)[grepl("sCo",rownames(a),fixed = T)]<-"stage"
  rownames(a)[grepl("aCol",rownames(a),fixed = T)]<-"age"
  rownames(a)[grepl("gCol",rownames(a),fixed = T)]<-"grade"
  a[,5]<-round(a[,5],digits = 4)
  colnames(a)<-c("lnhr","hr","se_lnhr","z","p")
  a<-as.data.frame(a)
  meta.i<-metagen(lnhr,se_lnhr,sm="HR",data = a,studlab = rownames(a))
  forest(meta.i, leftcols = c("studlab"),leftlabs = c("Study Element"),
         rightcols = c("effect.ci","p"),fontsize = 13,
         rightlabs = "p",just.studlab = "center",
         overall = F,col.square = "#363886",lwd=2,
         weight.study = "same",col.inside="black",xlab=names(selfEvalResult)[i])
}

HR<-sapply(selfEvalResult,function(x){x[[2]][1,2]})
lnHR<-sapply(selfEvalResult,function(x){x[[2]][1,1]})
selnHR<-sapply(selfEvalResult,function(x){x[[2]][1,3]})
low.HR<-exp(lnHR-1.96*selnHR)
high.HR<-exp(lnHR+1.96*selnHR)
C_index<-sapply(selfEvalResult,function(x){x[[3]][3]})
se_C_index<-sapply(selfEvalResult,function(x){x[[3]][4]})
low_c_index<-C_index-1.96*se_C_index
high_c_index<-C_index+1.96*se_C_index
HR_C<-data.frame(HR,lnHR,selnHR,low.HR,high.HR,
                 C_index,se_C_index,low_c_index,high_c_index)
lnHR<-meta.summaries(HR_C[,2],HR_C[,3],method="fixed",logscale=T)
if(lnHR$het[3]>0.05){
  summary_HR<-c(lnHR$summary,lnHR$se.summary)
}else{
  lnHR<-meta.summaries(HR_C[,2],HR_C[,3],method="random",logscale=T)
  summary_HR<-c(lnHR$summary,lnHR$se.summary)
}

C<-meta.summaries(HR_C[,6],HR_C[,7],method="fixed",logscale=F)
if(C$het[3]>0.05){
  summary_C<-c(C$summary,C$se.summary)
}else{
  C<-meta.summaries(HR_C[,6],HR_C[,7],method="random",logscale=F)
  summary_C<-c(C$summary,C$se.summary)
}
ssgsea_colla_Test<-c(summary_HR,summary_C)

ol<-colorRampPalette(c("#363886","#9ECAD4"))


forest(metagen(HR_C$lnHR,HR_C$selnHR,sm="HR"),
       col.square =  ol(nrow(HR_C)),col.diamond = "#70C1B3",
       studlab=rownames(HR_C))
       


forest(metagen(HR_C$C_index,HR_C$se_C_index,sm="C"),
       xlim = c(0,1),col.square =  ol(nrow(HR_C)),col.diamond = "#70C1B3",
       studlab=rownames(HR_C),ref = 0.5)
       
clinicalTrain<-clinical[trainPo]
clinicalTest<-clinical[-trainPo]
clinicalMatrixTrain<-getSurvivalMatrix(clinicalTrain,"global","Patch et al.")
clinicalMatrixTest<-getSurvivalMatrix(clinicalTest,"global","Patch et al.")
par(mfrow=c(3,4))#name:globalROC
globalTrain<-globalEval(scale(score),"global Training",clinicalMatrixTrain,1,2)
globalTest<-globalEval(scale(score.test),"global Test",clinicalMatrixTest,1,2)
global<-rbind(globalTrain,globalTest)

clinicalTrain<-clinical[trainPo]
clinicalTest<-clinical[testPo]
clinicalIndepend<-clinical[independPo]
clinicalMatrixTrain<-getSurvivalMatrix(clinicalTrain,"global","Patch et al.")
clinicalMatrixTest<-getSurvivalMatrix(clinicalTest,"global","Patch et al.")

x1<-clinical[[independPo[1]]]
phe1<-cbind(x1$days_to_death,x1$vital_status)
x2<-clinical[[independPo[[2]]]]
phe2<-cbind(x2$donor_survival_time,x2$donor_vital_status)
phe2[,2]<-ifelse(phe2[,2]==2,1,0)
phe<-rbind(phe1,phe2)
rownames(phe)<-c(rownames(x1),rownames(x2))
colnames(phe)<-c("time","status")
phe[,2]<-ifelse(phe[,2]=="0","Alive","Dead")
clinicalMatrixIndepend<-phe


globalTrain<-globalEval(scale(score),"global Trainingsst",clinicalMatrixTrain,1,2)
globalTest<-globalEval(scale(score.test),"global Testsst",clinicalMatrixTest,1,2)
globalIndepend<-globalEval(scale(score.test.independ),
                           "global independsst",clinicalMatrixIndepend,1,2)
                           
clinicalAll<-rbind(clinicalMatrixTrain,clinicalMatrixTest,clinicalMatrixIndepend)
score.all<-c(score,score.test,score.test.independ)

globalall<-globalEval(scale(score.all),"global all samplesst",clinicalAll,1,2)
data<-as.data.frame(global[,7:11])
data1<-data.frame(value=unlist(data),
                  survivalTime=rep(1:5,each=2),
                  class=rep(rownames(data),times=5))
ggplot(data1,aes(survivalTime,value,group=class))+
  geom_line(size=1,aes(color=class))+
  geom_point(size=2.5,aes(color=class))+
  geom_hline(yintercept = 0.5,linetype="dashed")+
  theme_bw()+
  scale_color_manual(values=c("#E07A9A","#805BA6"),
                     labels=c("Test set","Training set"))+
  labs(x="Survival time (Years)",y="AUC",title="Time depedent ROC curve")+
  theme(panel.border = element_blank(),panel.grid=element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(color = "black",size=15),
        axis.text.y=element_text(color="black",size=15),
        plot.title = element_text(colour = "black", 
                                  size = 15, vjust = 1,hjust = 0.5),
        axis.title = element_text(size=15),
        legend.title = element_blank(),
        legend.justification=c(0,0),
        legend.position=c(0.01,0.05),
        legend.text=element_text(size=12,face="bold"))
global<-as.data.frame(global)
global_Value_se<-data.frame(global$lnHR,global$selnHR,
                            global$C_index,global$se_C_index)
                            
global1<-data.frame(rownames(global),
                    format(global$HR,digits=4),
                    paste(format(global$low.HR,digits=4),
                          format(global$high.HR,digits=4),sep="-"),
                    format(global$p.HR,digits=3,scientific = T),
                    format(global$threeYearSurvival,digits=4),
                    format(global$fiveYearSurvival,digits=4),
                    format(global$C_index,digits=4),
                    paste(format(global$low_C_index,digits=4),
                          format(global$high_C_index,digits=4),sep="-"))
names(global1)<-c("Variables","HR","95%CI of HR","P value of HR",
                  "three-Year Survival",
                  "five-Year Survival","C-index","95%CI of C-index")

global$p.HR<-format(global$p.HR,digits = 3,scientific = T)
meta.i<-metagen(lnHR,selnHR,sm="HR",data = global,studlab = rownames(global))
forest(meta.i, leftcols = c("studlab"),leftlabs = c("Study Element"),
       rightcols = c("effect.ci","p.HR"),fontsize = 13,just.studlab = "center",
       overall = F,col.square = "#363886",lwd=2,rightlabs = "P value",
       weight.study = "same",col.inside="black")  

meta.i<-metagen(C_index,se_C_index,sm="C-index",data = global,
                studlab = rownames(global))
forest(meta.i, leftcols = c("studlab"),leftlabs = c("Study Element"),
       rightcols = c("effect.ci"),fontsize = 13,just.studlab = "center",
       overall = F,col.square = "#363886",lwd=2,xlim = c(0,1),
       weight.study = "same",col.inside="black",ref = 0.5)  