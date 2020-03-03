plotKM <-
  function(y, strata,
           xlab=paste("Time (",.guess.time.unit(y), ")", sep=""),
           ylab="Survival (%)",main=NULL,censor.at=NULL,show.n.risk=TRUE,
           n.risk.step=ifelse(show.n.risk,
                              round(ifelse(is.null(censor.at),max(y[,1],
                                                                  na.rm=TRUE),censor.at)/60)*10,
                              NULL),
           pval.defined=NULL,chisq.dist=NULL,
           show.legend=TRUE, show.n = !show.n.risk,  
           legendpos="bottomleft",npos="bottomleft",
           HRpos=ifelse(is.null(main),"margintop","topright"),
           show.HR=ifelse(length(levels(strata))==2,TRUE,FALSE),
           show.PV=TRUE,inverse.HR=TRUE,col=NULL,cex.base=1,
           cex.lab=cex.base,cex.legend=cex.base, cex.HR=cex.base,
           cex.n.risk=cex.base*0.85,cex.n=cex.base,lwd=2,
           lty=1:length(levels(strata)),...) {
    if (!is.factor(strata)) stop("strata not a factor")
    if (!is.null(censor.at)) { 
      y[y[,1] > censor.at,2] = 0
      y[y[,1] > censor.at,1] = censor.at
    }    
    fit =survfit(y~strata)
    
    # For tertile stratification, give HR for High vs. low risk
    if (length(grep("Mid Risk", strata))>0)
      fdiff=survdiff(y~strata,subset=strata!="Mid Risk")
    else fdiff =survdiff(y~strata)
    
    if (!is.null(pval.defined)) pv=pval.defined
    else {
      # If empirical chisqr distribution is given, use that one for p-value 
      # estimation
      if (is.null(chisq.dist)) {
        pv = 1-pchisq(fdiff$chisq,length(fdiff$n)-1)
      } else {
        pv = sum(chisq.dist>fdiff$chisq)/length(chisq.dist)
        pv = max(pv,1/length(chisq.dist))
      }
    }
    
    # Use some nicer colors 
    if (is.null(col)) col = .niceColorsF(strata)
    
    # set the margins if numbers at risk are to be shown
    ng <- length(levels(strata))
    old.mar <- par("mar")
    on.exit(par(mar = old.mar))
    .xaxt = "s"
    .xlab = xlab
    if (show.n.risk) {
      par(mar = old.mar + c(ng, 3*cex.n.risk, 0, 0))
      .xaxt = "n"
      .xlab = ""
    }
    
    plot(fit, col=col,xlab=.xlab,xaxt=.xaxt,las=1,
         ylab=ylab, yscale=100,cex.lab=cex.lab,lwd=lwd,lty=lty,main=main,...)
    legend.labels <- levels(as.factor(strata))
    if (show.n) legend.labels <- paste(legend.labels, " (n = ", fit$n, ")",
                                       sep="")
    if (show.legend) legend(legendpos, legend.labels,lwd=lwd,lty=lty,
                            col=col,cex=cex.legend)
    else if (show.n) {
      ncoords <- .corner.coord(npos,"npos")
      l <- "n = "
      tmp <- .corner.label(label=l, x=ncoords[1], y = ncoords[2],cex=cex.n) 
      for (i in 1:ng) {
        text(tmp$x+strwidth(l,cex=cex.n),tmp$y,adj=tmp$adj, label=fit$n[i],
             col=col[i],cex=cex.n)
        l <- paste(l, fit$n[i])
        if (i<ng)text(tmp$x+strwidth(l,cex=cex.n),tmp$y,adj=tmp$adj, 
                      label="/ ", cex=cex.n)
        l <- paste(l, "/")
      }
    }
    
    
    hr = ""
    if (show.HR || show.PV) {
      if (show.HR) {
        conf =summary(coxph(y~strata))$conf.int[1,]
        if (inverse.HR) { 
          conf = 1/conf
          conf[c(3,4)] = conf[c(4,3)]
        } 
        conf = sapply(conf, function(x) format(x, digits=3))
        hr = paste("HR ",conf[1],"; 95% CI, ", conf[3], " to ", 
                   conf[4], sep="")
      }
      if (show.HR & show.PV) hr = paste(hr, ";", sep="")
      if (show.PV) hr = paste(hr,"P",.format.pval.x(pv, digits=3));
      if (is.null(HRpos) || HRpos=="margintop") mtext(hr,side=3,cex=cex.HR)
      else {
        if (length(HRpos) > 1) { 
          text(HRpos[1], HRpos[2], hr, cex=cex.HR)
        }
        else { 
          HRcoords <- .corner.coord(HRpos,"HRpos")
          .corner.label(label=hr, x=HRcoords[1], y = HRcoords[2],
                        cex=cex.HR)        
        }
      }
    }
    
    if (show.n.risk) {
      require(survcomp)
      usr.xy <- par("usr")
      nrisk <- no.at.risk(Surv(surv,deceased)~strata,
                          data.frame(surv=y[,1],deceased=y[,2],strata), rep(TRUE,
                                                                            length(strata)), n.risk.step, floor(usr.xy[2]))
      at.loc <- seq(0, usr.xy[2], n.risk.step)
      axis(1, at = at.loc)
      mtext(xlab, side = 1, line = 2,cex=cex.lab)
      mtext("No. At Risk", side = 1, line = 3, at = -0.5 * 
              n.risk.step, adj = 1, cex = cex.n.risk, font = 2)
      for (i in 1:nrow(nrisk)) {
        mtext(levels(strata)[i], side = 1, line = 3 + i, at = -0.5 * 
                n.risk.step, adj = 1, cex = cex.n.risk)
        mtext(nrisk[i, -1], side = 1, at = at.loc, line = 3 + 
                i, adj = 1, cex = cex.n.risk)
      }
    }
    invisible(hr)
  }

.defaultLabels <- function(cutpoints) {
  if(length(cutpoints)==2) return( c("High Risk", "Mid Risk", "Low Risk"))
  c("High Risk", "Low Risk")
}

.guess.time.unit <- function(y) {
  m  = median(y[,1], na.rm=TRUE)
  mx = max(y[,1], na.rm=TRUE)
  
  if (m < 4 && mx < 25) return("Years")
  if (m > 200 || mx > 500) return("Days")
  "Months"
}

.format.pval.x <- function(pv, ...) {
  ifelse(pv < 0.001, "< 0.001", paste("=", format.pval(pv, ...)))
}

.corner.coord <- function(x, param.name) {
  switch(x,
         "topright"    = c(1,1),
         "topleft"     = c(-1,1),
         "bottomright" = c(1,-1),
         "bottomleft"  = c(-1,-1),
         stop(paste("Value of",param.name,"not supported")))
}

# taken from the plotrix package
.corner.label <- function(label = NULL, x = -1, y = 1,...)
{
  xoff <- strwidth("m")/2
  yoff <- strheight("m")/2
  par.usr <- par("usr")
  xpos <- par.usr[(3 + x)/2]
  ypos <- par.usr[(3 + y)/2 + 2]
  adj <- c((1 +  x)/2, (1 + y)/2)
  xpos <- xpos - x * xoff
  ypos <- ypos - y * yoff
  if (!is.null(label)) {
    text(xpos, ypos, label, adj = adj, ...)
  }
  return(list(x = xpos, y = ypos, adj=adj))
}


.niceColorsF <- function(strata,
                         palette=ifelse(is.factor(strata),"Set1","Greens"), 
                         n=ifelse(is.factor(strata),length(levels(strata)),8)) {
  max.n = switch(palette,
                 "Paired"  = 12,
                 "Pastel1" =  9,
                 "Set1"    =  9,
                 "Set3"    = 12, 8) 
  require(RColorBrewer)    
  if (n>max.n) nicecolors = colorRampPalette(brewer.pal(max.n, palette))(n)
  else nicecolors = brewer.pal(max(3,n),palette)    
  nicecolors
}
