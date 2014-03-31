\name{qqplot_JH}
\alias{qqplot_JH}
\title{Plots the pvalue quantiles against the uniform distribution. If adjusted pvalues are present, colors those that pass cutoff red}
\usage{
qqplot_JH(pvaldf, pval_colid = "pval", adj_pval_colid = "padj",
  adj_pval_cutoff = 0.05)
}
\arguments{
  \item{pvaldf}{dataframe, required}

  \item{pval_colid}{character string, column header of
  column containing unadjusted pvalues, required}

  \item{adj_pval_colid}{character string, column header of
  column containing adjusted pvalues, required}

  \item{adj_pval_cutoff,}{optional, defaults to 0.05}
}
\description{
Plots the pvalue quantiles against the uniform
distribution. If adjusted pvalues are present, colors those
that pass cutoff red
}
\examples{
qqplot_JH(df, pval_colid="pval", adj_pval_colid="padj", adj_pval_cutoff=0.05)
}

