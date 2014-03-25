\name{annotate_df}
\alias{annotate_df}
\title{Annotate Dataframe containing Ensembl IDs with Gene Symbols and Descriptions}
\usage{
annotate_df(df, df_ensemblid_header, biomart_ensembl_dataset,
  biomart_ensemblid_filter, biomart_genesymbol_attribute)
}
\arguments{
  \item{df}{dataframe, required}

  \item{df_ensemblid_header,}{character string containing
  the name of the dataframe column with ensembl ids,
  required}

  \item{biomart_ensembl_dataset,}{character string
  describing the biomart dataset to use, required}

  \item{biomart_ensemblid_filter,}{character string
  describing the biomart ensembl id filter to use,
  required}

  \item{biomart_genesymbol_attribute,}{character string
  describing the biomart gene symbol attribute to obtain,
  required}
}
\value{
annotated dataframe
}
\description{
Annotate Dataframe containing Ensembl IDs with Gene Symbols
and Descriptions
}
\examples{
annotate_df(temp, "id", 'mmusculus_gene_ensembl', "ensembl_gene_id", "mgi_symbol")
}
\seealso{
\code{\link{biomaRt}} used to annotate the dataframe
}

