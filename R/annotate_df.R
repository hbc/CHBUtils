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
  ensembl = useMart('ENSEMBL_MART_ENSEMBL', dataset = biomart_ensembl_dataset, host="www.ensembl.org")
  annot.df = getBM(attributes=c(biomart_ensemblid_filter, biomart_genesymbol_attribute, "description"),
            filters=c(biomart_ensemblid_filter), values=df[, df_ensemblid_header],
            mart=ensembl)
  m = merge(df, annot.df, by.x=df_ensemblid_header, by.y=biomart_ensemblid_filter)
  return(m)
}


get_biomart_transcript = function(v, type="description")
{
  require(biomaRt)
  mart = useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", host="www.ensembl.org")
  g <- getBM( attributes=c("ensembl_gene_id",type) , filters=
                "ensembl_gene_id"    , values =as.character(v) ,mart=mart)
  #row.names(g) = g[,1]
  out <- c("", length=length(v))
  idx <- match(g[,1], v)
  out[idx] <- g[,2]
  names(out[idx]) <- g[,1]
  out
  #g
}


# convertIDs(row.names(.res), "ENSEMBLTRANS", "SYMBOL", org.Mm.eg.db, "useFirst")
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

#rewrite plotsplice function
plotSplice =
  function (fit, top, coef = ncol(fit),  geneid = NULL, genecolname = NULL,
            rank = 1L, FDR = 0.05)
  {
    if (is.null(genecolname))
      genecolname <- fit$genecolname
    else genecolname <- as.character(genecolname)
    if (is.null(geneid)) {
      if (rank == 1L)
        i <- which.min(fit$gene.F.p.value[, coef])
      else i <- order(fit$gene.F.p.value[, coef])[rank]
      geneid <- paste(fit$gene.genes[i, genecolname], collapse = ".")
    }
    else {
      geneid <- as.character(geneid)
      i <- which(fit$gene.genes[, genecolname] == geneid)[1]
      if (!length(i))
        stop(paste("geneid", geneid, "not found"))
    }
    j <- fit$gene.firstexon[i]:fit$gene.lastexon[i]
    exoncolname <- fit$exoncolname
    strcol <- grepl("strand", colnames(fit$gene.genes), ignore.case = TRUE)
    if (any(strcol))
      geneid <- paste0(geneid, " (", as.character(fit$gene.genes[i,
                                                                 strcol])[1], ")")
    if (is.null(exoncolname)) {
      plot(fit$coefficients[j, coef], xlab = "Exon", ylab = "logFC (this exon vs rest)",
           main = geneid, type = "b")
    }
    else {
      exon.id <- fit$genes[j, exoncolname]
      xlab <- paste("Exon", exoncolname, sep = " ")
      plot(fit$coefficients[j, coef], xlab = "", ylab = "logFC (this exon vs rest)",
           main = geneid, type = "b", xaxt = "n", cex=1)
      axis(1, at = 1:length(j), labels = exon.id, las = 2,
           cex.axis = 1)
      # mtext(xlab, side = 1, padj = 5.2)
      # top <- topSplice(fit, coef = coef, number = Inf, test = "t",
      #                  FDR = FDR)

      m <- which(top[, genecolname] %in% as.character(fit$gene.genes[i,genecolname]))
      mark = NULL
      if (length(m) > 0) {
        if (length(m) == 1)
          cex <- 1.5
        else {
          abs.fdr <- abs(log10(top$FDR[m]))
          from <- range(abs.fdr)
          to <- c(1, 2)
          cex <- (abs.fdr - from[1])/diff(from) * diff(to) +
            to[1]
        }
        mark <- match(top[m, exoncolname], exon.id)
        points((1:length(j))[mark], fit$coefficients[j[mark],
                                                     coef], col = "red", pch = 16, cex = cex)
      }
    }
    abline(h = 0, lty = 2)
    invisible()
    mark
  }

# reduce clusterprofile output
reduce_cp = function(genes, lim=100){
  seen = c()
  idx = sapply(genes, function(x){
    if (x == "")
      return(FALSE)
    here = as.character(as.vector(unlist(strsplit(x, split = "/"))))
    if (length(here) > lim)
      return(FALSE)
    c = intersect(seen, here)
    seen <<- unique(c(seen, here))
    score = 0.9 * length(here)
    if (length(c) < score)
      return(TRUE)
    return(FALSE)
  })
  idx
}

.isvalid= function(dd){
  if (is.null(dim(dd))){
    return(FALSE)
  }
  if (nrow(dd)==0){
    return(FALSE)
  }
  return(TRUE)
}

