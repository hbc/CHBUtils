\name{mdsplot.RGset}
\alias{mdsplot.RGset}
\title{MDS (MultiDimensional Scaling) Plots}
\usage{
mdsplot.RGset(dat, numPositions = 1000, sampGroups, sampNames,
  main = "MDS plot", alpha = 1, cex = 1)
}
\arguments{
  \item{dat}{object, RGChannelSet, required}

  \item{numPositions,}{integer, number of most variable
  positions to use for MDS, defaults to 1000}

  \item{sampGroups,}{chracter string, column name in
  phenoData of RGset containing groups for coloring groups
  on plot, required}

  \item{sampNames,}{chracter string , column name in
  phenoData of RGset containing sample labels for plot,
  optional}

  \item{main}{chracter string, title for the plot, defaults
  to "MDSplot"}

  \item{alpha}{floating point, alpha shading from 0
  (transparent) to 1 (opaque), defaults to 1}

  \item{cex}{integer, size of points in plot, defaults to
  1}
}
\value{
plot of samples plotted by MDS
}
\description{
MDS (MultiDimensional Scaling) Plots
}
\examples{
mdsplot.RGset(dat=RGset, numPositions=1000, sampGroups="groups", sampNames="sampleIDs", main="MDSplot", alpha=0.8, cex=10)
}
\seealso{
\code{\link{cmdscale}} which this function uses to
calculate MDS
}

