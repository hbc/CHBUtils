\name{PCAplot.eset}
\alias{PCAplot.eset}
\title{Pairwise Principal Component Analysis Plots}
\usage{
PCAplot.eset(eset = NULL, categories = NULL, title = NULL,
  colorpalette = NULL, alpha = 1, numcomponents = 4)
}
\arguments{
  \item{eset}{ExpressionSet, required}

  \item{categories}{character vector column identifier
  containing different categories to highlight with colors,
  required}

  \item{secondary_categories}{character vector column
  identifier containing different categories to highlight
  by shape, optional}

  \item{title}{title for the plot, required}

  \item{colorpalette}{vector of RGB colors (long enough for
  a color per category), required}

  \item{alpha}{alpha shading from 0 (transparent) to 1
  (opaque), defaults to 1}

  \item{numcomponents}{number of principal components to
  plot, defaults to 4}
}
\value{
pairwise plots of samples plotted by principal component
}
\description{
Pairwise Principal Component Analysis Plots
}
\examples{
PCAplot.eset(eset=AffyNorm, categories="groups", title="PCAplot - groups", colorpalette=c("#FF0000", "#00FF00", "#0000FF", alpha=0.8, numcomponents=4)
}
\seealso{
\code{\link{pair}} which this function uses to plot all
pairwise combos
}

