#' Annotate Dataframe containing Ensembl IDs with Gene Symbols and Descriptions
#'
#' @param df dataframe, required
#' @param df_ensemblid_header, character string containing the name of the dataframe column with ensembl ids, required
#' @param biomart_ensembl_dataset, character string describing the biomart dataset to use, required
#' @param biomart_ensemblid_filter, character string describing the biomart ensembl id filter to use, required 
#' @param biomart_genesymbol_attribute, character string describing the biomart gene symbol attribute to obtain, required
#' @return annotated dataframe
#' @seealso \code{\link{biomaRt}} used to annotate the dataframe
#' @export
#' @examples
#' annotate_df(temp, "id", 'mmusculus_gene_ensembl', "ensembl_gene_id", "mgi_symbol")

annotate_df = function(df, df_ensemblid_header, biomart_ensembl_dataset, biomart_ensemblid_filter, biomart_genesymbol_attribute) {
  require(biomaRt)
  ensembl = useMart('ensembl', dataset = biomart_ensembl_dataset)
  annot.df = getBM(attributes=c(biomart_ensemblid_filter, biomart_genesymbol_attribute, "description"),
            filters=c(biomart_ensemblid_filter), values=df[, df_ensemblid_header],
            mart=ensembl)
  m = merge(df, annot.df, by.x=df_ensemblid_header, by.y=biomart_ensemblid_filter)
  return(m)
}
