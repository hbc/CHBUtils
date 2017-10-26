#' MDS (MultiDimensional Scaling) Plots
#'
#' @param dat object, RGChannelSet, required
#' @param numPositions, integer, number of most variable positions to use for MDS, defaults to 1000
#' @param sampGroups, chracter string, column name in phenoData of RGset containing groups for coloring groups on  plot, required
#' @param sampNames, chracter string , column name in phenoData of RGset containing sample labels for plot, optional
#' @param main chracter string, title for the plot, defaults to "MDSplot"
#' @param alpha floating point, alpha shading from 0 (transparent) to 1 (opaque), defaults to 1
#' @param cex integer, size of points in plot, defaults to 1 
#' @return plot of samples plotted by MDS
#' @seealso \code{\link{cmdscale}} which this function uses to calculate MDS
#' @export
#' @examples
#' mdsplot.RGset(dat=RGset, numPositions=1000, sampGroups="groups", sampNames="sampleIDs", main="MDSplot", alpha=0.8, cex=10)

mdsplot.RGset <- function (dat, numPositions = 1000, sampGroups, sampNames, main = "MDS plot", alpha=1, cex=1) {
  if (missing(dat)){
    stop("data required")
  }
  if (missing(sampGroups)){
    stop("sampGroups required")
  }
  # get beta values
  b <- getBeta(dat)
  pd <- pData(dat)
  # reorder and subset to most variable rows
  o <- order(-rowVars(b))[1:numPositions]
  # calculate euclidean distances  
  d <- dist(t(b[o, ]))
  # mds
  fit <- cmdscale(d)
  # prep for ggplot
  fit <- as.data.frame(fit)
  names(fit) <- c("x", "y")
  # merge in phenoDAta to mds fit data
  fit <- cbind(fit, pd)
  # generate color factors for ggplot
  fit$col <- fit[,sampGroups]
  fit$col <- as.character(fit$col)
  fit <- as.data.frame(fit) # Required again as new minfi returns DataFrame class, not data.frame
  # factor colors in palette by sampGroups
  if (missing(sampNames)) {
    ggplot(fit, aes(x=x, y=y, color=col))+
      geom_point(size=cex, alpha=alpha)+
      ggtitle(main)+
      scale_color_discrete(name="Experimental\nCondition")+
      scale_x_continuous(name="")+
      scale_y_continuous(name="")
  }  else {
    fit$labels <- factor(pd[,sampNames])
    ggplot(fit, aes(x=x, y=y, color=col))+
      geom_point(size=cex, alpha=alpha)+
      ggtitle(main)+
      geom_text(data = fit, aes(x,y, label = labels), hjust = 2)+
      scale_color_discrete(name="Experimental\nCondition")+
      scale_x_continuous(name="")+
      scale_y_continuous(name="")
  }
}
