#' Plots the differntially expressed genes on an MA plot, ggplot style
#' If adjusted pvalues are present, colors those that pass cutoff red
#' @param pvaldf dataframe, required
#' @param basemean_colid character string, column header of column containing baseMean values, required
#' @param log2foldchange_colid character string, column header of column containing log2 transformed fold change values, required
#' @param adj_pval_colid character string,column header of column containing adjusted pvalues, required
#' @param adj_pval_cutoff, optional, defaults to 0.05
#' @param title character string, title for the plot, optional, defaults to "MA plot"
#' @export
#' @examples
#' plotDE(df, basemean_colid="baseMean",log2foldchange_colid=log2FoldChange, adj_pval_colid="padj", adj_pval_cutoff=0.05, title="My MA plot")

plotDE = function(pvaldf,  basemean_colid=NULL, log2foldchange_colid=NULL, adj_pval_colid=NULL, adj_pval_cutoff=0.05, title="MA plot")  {
  title = paste("M-A plot of", sep = " - ")
  pvaldf$colors <- ifelse(pvaldf[,adj_pval_colid] < adj_pval_cutoff, "sig", "nonsig")
  pvaldf$bm <- pvaldf[,basemean_colid]
  pvaldf$lfc <- pvaldf[,log2foldchange_colid]
  plot <- ggplot(data = pvaldf, aes(x = log(bm), y = lfc, colour = colors)) + 
    geom_point(size = 3)  + 
    scale_colour_manual(name = "BFH adjusted pvalue", values = c("#00000033", "#FF0000FF"), labels = c(paste("q>", adj_pval_cutoff, sep = ""), paste("q<", adj_pval_cutoff, sep = ""))) + 
    labs(title = title)
  plot
}
