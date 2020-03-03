plotSurvivalCurve<-function(fit,data,title){
  ggsurvplot(fit, data =data, title=title,
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