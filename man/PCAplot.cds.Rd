\name{PCAplot.cds}
\alias{PCAplot.cds}
\title{Pairwise Principal Component Analysis Plots of CountDataSet objects}
\usage{
PCAplot.cds(countdataset, categories, title, colorpalette, alpha = 1,
  numcomponents = 4, normalize.counts = TRUE)
}
\arguments{
  \item{countdataset}{CountDataSet, required}

  \item{categories}{character string, column header of
  column containing different categories to highlight with
  colors, required}

  \item{title}{character string, title for the plot,
  required}

  \item{colorpalette}{charcter vector of RGB colors (long
  enough for a color per category), required}

  \item{alpha}{alpha shading from 0 (transparent) to 1
  (opaque), defaults to 1}

  \item{numcomponents}{number of principal components to
  plot, defaults to 4}
}
\value{
pairwise plots of samples plotted by principal component
}
\description{
Pairwise Principal Component Analysis Plots of CountDataSet
objects
}
\examples{
PCAplot.cds(countdataset=cds, categories="condition", title="PCAplot", colorpalette=c("#FF0000", "#00FF00", "#0000FF"), alpha=0.8, numcomponents=4)
}
\seealso{
\code{\link{pair}} which this function uses to plot all
pairwise combos
}

