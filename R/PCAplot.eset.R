#' Pairwise Principal Component Analysis Plots
#'
#' @param eset ExpressionSet, required
#' @param categories character vector column identifier containing different categories to highlight with colors, required
#' @param secondary_categories character vector column identifier containing different categories to highlight by shape, optional
#' @param title title for the plot, required
#' @param colorpalette vector of RGB colors (long enough for a color per category), required
#' @param alpha alpha shading from 0 (transparent) to 1 (opaque), defaults to 1
#' @param numcomponents number of principal components to plot, defaults to 4 
#' @return pairwise plots of samples plotted by principal component
#' @seealso \code{\link{pair}} which this function uses to plot all pairwise combos
#' @export
#' @examples
#' PCAplot.eset(eset=AffyNorm, categories="groups", title="PCAplot - groups", colorpalette=c("#FF0000", "#00FF00", "#0000FF", alpha=0.8, numcomponents=4)
  

PCAplot.eset <- function(eset=NULL, categories=NULL, secondary_categories=NULL, title=NULL, colorpalette=NULL, alpha=1, numcomponents=4){
  if (missing(eset)){
    stop("Function requires an ExpressionSet")
  }
  if(missing(categories)){
    stop("Function requires you to specify phenoData column name (as categories to highlight)")
  }
  if(missing(title)){
    stop("Function requires you to specify a title for plot")
  }
  if (missing(colorpalette)){
    stop("Funciton requires you to specify a colorpalette")
  }
  if (missing(secondary_categories)){
    shapes=21
  } else {
    shapes <- c(21:25)[factor(as.character(unlist(pd[,secondary_categories])))]
  }
  alpha <- sprintf("%x", ceiling(alpha*255))
  colorpalette <- paste(colorpalette, alpha, sep="")
  eset.core <- exprs(eset) 
  myPca.core <- prcomp(t(eset.core))
  tmpPCAData.core <- as.data.frame(myPca.core$x[,1:numcomponents])
  pd <- pData(eset)
  colors <- colorpalette[factor(as.character(unlist(pd[,categories])))]
  if(missing(secondary_categories)){
    legend_values <- unique(cbind(colors, as.character(pd[,categories])))
  } else {
    legend_values <- unique(do.call(cbind, list(shapes, colors, as.character(pd[,categories]), as.character(pd[,secondary_categories]))))
  }
  
  legend_values=unique(cbind(colors, as.character(pd[,categories])))
  pairs(tmpPCAData.core, bg=colors, col="#606060", cex=2, pch=21, main=title, oma=c(8,5,5,14))
  legend("right", cex=0.7, col="#606060", pt.bg=legend_values[,1], pt.cex=1.5, legend=legend_values[,2],  pch=21, bty="n", x.intersp=1)
}
