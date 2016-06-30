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
  ensembl = useMart('ENSEMBL_MART_ENSEMBL', dataset = biomart_ensembl_dataset, host="useast.ensembl.org")
  annot.df = getBM(attributes=c(biomart_ensemblid_filter, biomart_genesymbol_attribute, "description"),
            filters=c(biomart_ensemblid_filter), values=as.character(df[, df_ensemblid_header]),
            mart=ensembl)
  m = merge(df, annot.df, by.x=df_ensemblid_header, by.y=biomart_ensemblid_filter, all.x=T)
  return(m)
}


#' Get biomart annotation from id
#'
#' @param v vector of ids
#' @param type character name to get from biomart
#' @param id character defining what kind of id is provided
#' @param dataset specie database to connect to biomart
#' @return annotated vector in the same order
#' @seealso \code{\link{biomaRt}} used to annotate the dataframe
#' @export
#' @examples
#' get_biomrt(v, "gene_biotype")
get_biomart = function(v, type="description", id="ensembl_gene_id",
                       dataset="hsapiens_gene_ensembl")
{
  require(biomaRt)
  v = as.character(v)
  v[is.na(v)] = "None"
  mart = useMart("ENSEMBL_MART_ENSEMBL", dataset, host="www.ensembl.org")
  g <- getBM( attributes=c(id,type) , filters=
                id, values =as.character(v) ,mart=mart)
  out <- rep("", length(v))
  idx <- match(g[,1], v)
  out[idx] <- g[,2]
  names(out[idx]) <- g[,1]
  out
}

#' Get other Ids for genes in the same order
#'
#' @param ids, character vector of gene ids
#' @param from, type of ID of ids vector
#' @param to, type of ID to convert
#' @param db, org.species.eg.db class
#' @param ifMultiple, what to do when multiple possibles IDs
#' @return annotated dataframe
#' @seealso \code{\link{AnnotationDbi::select}} used to annotate the dataframe
#' @export
#' @examples
#' convertIDs(list_genes, "ENSEMBLTRANS", "SYMBOL", org.Mm.eg.db, "useFirst")
convertIDs <- function( ids, from, to, db, ifMultiple=c("putNA", "useFirst")) {
  stopifnot( inherits( db, "AnnotationDb" ) )
  ifMultiple <- match.arg( ifMultiple )
  if (sum(ids %in% keys(db, from))==0)
    return(ids)
  suppressMessages( selRes <- AnnotationDbi::select(
    db, keys=ids, keytype=from, columns=c(from,to) ) )
  if ( ifMultiple == "putNA" ) {
    duplicatedIds <- selRes[ duplicated( selRes[,from] ), from ]
    selRes <- selRes[ ! selRes[,from] %in% duplicatedIds, ]
  }
  return( selRes[ match( ids, selRes[,from] ), to ] )
}


