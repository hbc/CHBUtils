% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/boxplotoutliers.R
\name{boxplotoutliers}
\alias{boxplotoutliers}
\title{Cel file IDs of intensity boxplot outlier arrays (in characters).}
\usage{
boxplotoutliers(eset = NULL, logtransform = TRUE)
}
\arguments{
\item{eset}{GeneFeatureSetm, AffyBatch or ExpressionSet object}

\item{logtransform}{Boolean, log transform data before processing?}
}
\value{
character vector of arrays (CEL file IDs) that are outliers for intensity boxplots.
}
\description{
Cel file IDs of intensity boxplot outlier arrays (in characters).
}
\examples{
boxplotoutliers(eset=AffyRaw, logtransform=TRUE)s
}
\seealso{
\code{\link{arrayQualityMetrics}} which this function subsets
}

