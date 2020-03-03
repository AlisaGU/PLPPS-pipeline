Bubble<-function(data,title,fillStr,gradientHigh,gradientLow,hueFirst,hueLast){
  #data contains X,Y,Similarity,P,Significance
  #names are separately D1,D2,Similarity,p,Significance
  data$Significance = factor(data$Significance)
  if(fillStr=="fillAndColour"){
    ggplot(data,aes(D1,D2))+
      geom_point(shape = 21, aes(fill=p,size=Similarity,colour=Significance))+
      scale_fill_continuous(high = gradientHigh,low = gradientLow) +
      scale_color_manual(values = c(hueFirst,hueLast))+
      theme_bw() +
      theme(panel.border = element_blank(),panel.grid=element_line(),
            axis.line = element_line(colour = "black"),
            axis.text.x = element_text(angle = 30,hjust = 1,vjust = 1),
            plot.title = element_text(colour = "black", face = "bold", 
                                      size = 14, vjust = 1,hjust = 0.5),
            axis.title.x = element_blank(),
            axis.title.y = element_blank())+
      labs(title=title)
  }
  if(fillStr=="colour"){
    ggplot(data,aes(D1,D2))+
      geom_point(shape = 16, aes(size=Similarity,colour=Significance))+
      scale_color_manual(values = c(hueFirst,hueLast))+
      theme_bw() +
      theme(panel.border = element_blank(),panel.grid=element_line(),
            axis.line = element_line(colour = "black"),
            axis.text.x = element_text(angle = 30,hjust = 1,vjust = 1),
            plot.title = element_text(colour = "black", face = "bold", 
                                      size = 14, vjust = 1,hjust = 0.5),
            axis.title.x = element_blank(),
            axis.title.y = element_blank())+
      labs(title=title)
  }
  if(fillStr=="fill"){
    ggplot(data,aes(D1,D2))+
      geom_point(shape = 21,aes(size=Similarity,fill=Similarity),color="white")+
      #scale_fill_continuous(low="#dadaeb",high="#6a51a3")+#PUPPLE
      #scale_fill_continuous(low="#fdd0a2",high="#d94801")+#cheng
      scale_fill_continuous(low="#c6dbef",high="#2171b5")+#blue
      #scale_fill_continuous(low="#fff5f0",high="#cb181d")+#red
      geom_point(shape=21,aes(colour=Significance,size=Similarity))+
      scale_color_manual(values = c("#FFFFFF","#FF0000"))+
      scale_size_continuous(guide = guide_legend(reverse=TRUE))+
      theme_bw() +
      theme(panel.border = element_blank(),panel.grid=element_blank(),
            axis.line = element_line(colour = "black"),
            axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1,color="black"),
            axis.text.y = element_text(color="black"),
            plot.title = element_text(colour = "black", face = "bold", 
                                      size = 14, vjust = 1,hjust = 0.5),
            axis.title.x = element_blank(),
            axis.title.y = element_blank())+
      labs(title=title)
  }
}