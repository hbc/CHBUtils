\name{plotDE}
\alias{plotDE}
\title{Plots the differntially expressed genes on an MA plot, ggplot style
If adjusted pvalues are present, colors those that pass cutoff red}
\usage{
plotDE(pvaldf, basemean_colid = NULL, log2foldchange_colid = NULL,
  adj_pval_colid = NULL, adj_pval_cutoff = 0.05)
}
\arguments{
  \item{pvaldf}{dataframe, required}

  \item{basemean_colid}{character string, column header of
  column containing baseMean values, required}

  \item{log2foldchange_colid}{character string, column
  header of column containing log2 transformed fold change
  values, required}

  \item{adj_pval_colid}{character string,column header of
  column containing adjusted pvalues, required}

  \item{adj_pval_cutoff,}{optional, defaults to 0.05}
}
\description{
Plots the differntially expressed genes on an MA plot,
ggplot style If adjusted pvalues are present, colors those
that pass cutoff red
}
\examples{
plotDE(df, basemean_colid="baseMean",log2foldchange_colid=log2FoldChange, adj_pval_colid="padj", adj_pval_cutoff=0.05)
}

