#' Cel file IDs of intensity boxplot outlier arrays (in characters).
#'
#' @param eset GeneFeatureSetm, AffyBatch or ExpressionSet object
#' @param logtransform Boolean, log transform data before processing?
#' @return character vector of arrays (CEL file IDs) that are outliers for intensity boxplots.
#' @seealso \code{\link{arrayQualityMetrics}} which this function subsets
#' @export
#' @examples
#' boxplotoutliers(eset=AffyRaw, logtransform=TRUE)s

boxplotoutliers <- function(eset=NULL, logtransform=TRUE){
  require(arrayQualityMetrics)
  preparedData=prepdata(expressionset=eset, intgroup=NULL, do.logtransform=TRUE)
  boxplot<-aqm.boxplot(preparedData)
  outliers <- names(boxplot@outliers@which)  
  return(outliers)
}
