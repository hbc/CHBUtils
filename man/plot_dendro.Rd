\name{plot_dendro}
\alias{plot_dendro}
\title{ggplot based dendrogram plot}
\usage{
plot_dendro(x, title = "", labels.colname = NULL, colors.colname = NULL)
}
\arguments{
  \item{x}{ExpressionSet, required}

  \item{title,}{character vector with title of the plot,
  optional}

  \item{labels.colname}{character vector of length one with
  column id of categorical labels, required}

  \item{colors.colname}{character vector of length one with
  column id of categorical colors, required}
}
\value{
dendrogram plot with categorical labels and colors
}
\description{
ggplot based dendrogram plot
}
\examples{
plot_dendro <- function(eset, title="mydendrogram", labels.colname="age", colors.colname="genotype")
}
\seealso{
\code{\link{ggdendro}} which this function modifies
}

