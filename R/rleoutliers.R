#' Cel file IDs of intensity RLE outlier arrays (in characters).
#'
#' @param eset AffyBatch
#' @param logtransform Boolean, log transform data before processing?
#' @return character vector of arrays (CEL file IDs) that are outliers for intensity boxplots.
#' @seealso \code{\link{arrayQualityMetrics}} which this function subsets
#' @export
#' @examples
#' rleoutliers(eset=AffyRaw, logtranform=TRUE)

rleoutliers <- function(eset=NULL, logtransform=TRUE){
  if(class(eset)!="AffyBatch"){
    stop("This function requires microarray data in the AffyBatch format (i.e. from 3'UTR arrays like the Affymetrix U133 or 430A series")
  } else {
    require(arrayQualityMetrics)
    preparedData=prepdata(eset, intgroup=NULL, do.logtransform=logtransform)
    preparedAffy=prepaffy(eset, preparedData)
    rle<-aqm.rle(preparedAffy)
    outliers <- names(nuse@outliers@which)  
  }
  return(outliers)
}