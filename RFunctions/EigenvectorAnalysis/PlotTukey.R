## Functions to plot the difference in mean values and confidence intervals + p value for tukey test 
## Lauren Overend
## lauren.overend@oriel.ox.ac.uk
## September 2022
## Adapted from: https://rpubs.com/brouwern/plotTukeyHSD2


#### The function STARTS here ####
plotTukeyHSD <- plotTukeysHSD <- function(tukey.out, modulename, day, y.min=NA, y.max=NA, x.axis.label = "Comparison", y.axis.label = "Difference in Mean", axis.adjust = 0, adjust.x.spacing = 5){
  tukey.out <- as.data.frame(tukey.out)
  means <- tukey.out$diff
  categories <- row.names(tukey.out)
  groups <- length(categories)
  ci.low <- tukey.out$lwr
  ci.up  <- tukey.out$upr                         
  n.means <- length(means)
  tukey.out$sig <- NA
  tukey.out[,"p adj"] <- as.numeric(tukey.out[,"p adj"])
  tukey.out$sig[tukey.out[,"p adj"]>=0.1] <- "ns"
  tukey.out$sig[tukey.out[,"p adj"]<0.1 & tukey.out[,"p adj"]>=0.05] <- "."
  tukey.out$sig[tukey.out[,"p adj"]<0.05 & tukey.out[,"p adj"]>=0.01] <- "*"
  tukey.out$sig[tukey.out[,"p adj"]<0.01 & tukey.out[,"p adj"]>=0.001] <- "**"
  tukey.out$sig[tukey.out[,"p adj"]<0.001  & tukey.out[,"p adj"]>=0.0001] <- "***"
  tukey.out$sig[tukey.out[,"p adj"]<0.0001] <- "****"
  significance <- tukey.out[,c("diff", "sig")]
  significance$x <- rownames(significance)

  #determine where to plot points along x-axis
  x.values <- 1:n.means
  x.values <- x.values/adjust.x.spacing                         
  # calculate values for plotting limits            
  
  if(is.na(y.min)){
	  y.max <- max(ci.up) +                    
		max(ci.up)*axis.adjust
	  y.min <- min(ci.low) - 
		max(ci.low)*axis.adjust
	  if(y.min>0){
		y.min <- 0
	 }
  } else {
	y.min=y.min
	y.max=y.max
  }
  if(groups == 2){ x.values <- c(0.25, 0.5)}
  if(groups == 3){ x.values <- c(0.25, 0.5,0.75)}
  x.axis.min <- min(x.values)-0.05
  x.axis.max <- max(x.values)+0.05
  x.limits <- c(x.axis.min,x.axis.max)
  significance$y <- (y.max+0.1)
  #Plot means
  plot(means ~ x.values,
       xlim = x.limits,
       ylim = c(y.min,(y.max+0.1)),
       xaxt = "n",
       xlab = "",
       ylab = "",
       cex = 1.25,
	   main=paste0(modulename, "\n", day),
       pch = 16)
  axis(side = 1, 
       at = x.values,
       labels = categories,
      )
  #Plot upper error bar 
  lwd. <- 2
  arrows(y0 = means,
         x0 = x.values,
         y1 = ci.up,
         x1 = x.values,
         length = 0,
         lwd = lwd.)
  #Plot lower error bar
  arrows(y0 = means,
         x0 = x.values,
         y1 = ci.low,
         x1 = x.values,
         length = 0,
         lwd = lwd.) 
  #add reference line at 0
  abline(h = 0, col = 2, lwd = 2, lty =2)
  mtext(text = x.axis.label,side = 1,line = 2.5)
  mtext(text = y.axis.label,side = 2,line = 2.5, cex=1)
  mtext(text = "Error bars = 95% CI",side = 3,line = 0,adj = 0, cex=2/3)
  text(x=x.values, y=significance$y, significance[,2], cex=2.5, col="red")
  return(c(y.min, y.max))
  
}