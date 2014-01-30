#' Cel file IDs of intensity NUSE outlier arrays (in characters).
#'
#' @param eset AffyBatch
#' @param logtransform Boolean, log transform data before processing?
#' @return character vector of arrays (CEL file IDs) that are outliers for intensity boxplots.
#' @seealso \code{\link{arrayQualityMetrics}} which this function subsets
#' @export
#' @examples
#' nuseoutliers(eset=AffyRaw, logtranform=TRUE)

nuseoutliers <- function(eset=NULL, logtransform=TRUE){
  if(class(eset)!="AffyBatch"){
    stop("This function requires microarray data in the AffyBatch format (i.e. from 3'UTR arrays like the Affymetrix U133 or 430A series")
  } else {
    require(arrayQualityMetrics)
    preparedData=prepdata(eset, intgroup=NULL, do.logtransform=logtransform)
    preparedAffy=prepaffy(eset, preparedData)
    nuse<-aqm.nuse(preparedAffy)
    outliers <- names(nuse@outliers@which)  
  }
  return(outliers)
}

