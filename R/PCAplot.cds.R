#' Pairwise Principal Component Analysis Plots of CountDataSet objects
#'
#' @param countdataset CountDataSet, required
#' @param categories character string, column header of column containing different categories to highlight with colors, required
#' @param title character string, title for the plot, required
#' @param colorpalette charcter vector of RGB colors (long enough for a color per category), required
#' @param alpha alpha shading from 0 (transparent) to 1 (opaque), defaults to 1
#' @param numcomponents number of principal components to plot, defaults to 4 
#' @return pairwise plots of samples plotted by principal component
#' @seealso \code{\link{pair}} which this function uses to plot all pairwise combos
#' @export
#' @examples
#' PCAplot.cds(countdataset=cds, categories="condition", title="PCAplot", colorpalette=c("#FF0000", "#00FF00", "#0000FF"), alpha=0.8, numcomponents=4)

PCAplot.cds <- function(countdataset, categories, title, colorpalette, alpha=1, numcomponents=4, normalize.counts=TRUE){
  # check for missing options
  if (missing(countdataset)){
    stop("Function requires an CountDataSet")
  }
  if(missing(categories)){
    stop("Function requires you to specify phenoData column name (as categories to highlight)")
  }
  if(missing(title)){
    stop("Function requires you to specify a title for plot")
  }
  if (missing(colorpalette)){
    stop("Function requires you to specify a colorpalette")
  }
  # get metadata
  pd <- pData(countdataset)
  # set colors and alphas
  alpha <- sprintf("%x", ceiling(alpha*255))
  colorpalette <- paste(colorpalette, alpha, sep="")
  colors <- colorpalette[factor(as.character(unlist(pd[,categories])))]
  # normalize data if requested
  if (normalize.counts){
    df <- counts(countdataset, normalized=TRUE)
  } else {
    df <- counts(countdataset, normalized=FALSE)
  }
  # run the actual PCA
  myPca.core <- prcomp(t(df))
  tmpPCAData.core <- as.data.frame(myPca.core$x[,1:numcomponents])
  # plot the result
  pairs(tmpPCAData.core, bg=colors, col="#606060", cex=2, pch=21, main=title, oma=c(8,5,5,14))
  legend("right", cex=0.7, col="#606060", pt.bg=colors, pt.cex=1.5, legend=pd[,categories],  pch=21, bty="n", x.intersp=1)
}
  
