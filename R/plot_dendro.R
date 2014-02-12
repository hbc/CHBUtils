#' ggplot based dendrogram plot 
#'
#' @param x ExpressionSet, required
#' @param title, character vector with title of the plot, optional
#' @param labels.colname character vector of length one with column id of categorical labels, required
#' @param colors.colname character vector of length one with column id of categorical colors, required
#' @return dendrogram plot with categorical labels and colors
#' @seealso \code{\link{ggdendro}} which this function modifies
#' @export
#' @examples
#' plot_dendro <- function(eset, title="mydendrogram", labels.colname="age", colors.colname="genotype")

plot_dendro <- function(x, title="", labels.colname=NULL, colors.colname=NULL) {
  require(ggdendro)
  meta.x <- pData(x)
  # force the metadata into character format so you don't end up with gradient/continuous color schemes for numerical variables in the final plot  
  meta.x <- as.matrix(meta.x) 
  ## do the actual statistics and put into dendrogram 
  myDist <- dist(t(exprs(x)))
  myTree <-hclust(myDist)
  dhc <- as.dendrogram(myTree)
  ddata <- dendro_data(dhc, type="rectangle")
  # the labels of the dendrogram are pulled from the Expression set exprs column names, it's nice to rename them to something more intelligible if you haven't already, as well as match them up to the metadata for label coloring
  ## check to see if the column names of the expression set match anything in the metadata, or match the rownames
  if (identical(colnames(exprs(x)), row.names(meta.x))) {
    meta.x <- row2colnames(meta.x, "rownames")
    matchcol <- "rownames"
  } else if (any(apply(meta.x, 2, function(column) identical(as.character(unlist(column)), colnames(exprs(x)))))) {
    matchcol <- names(which(apply(meta.x, 2, function(column) identical(as.character(unlist(column)), colnames(exprs(x))))))
  } else {
    print("ExpressionSet sampleNames and pData row.names or pData column must match")
    stop()
  }
  ## merge the metadata with the dendrogram labels using the commmon column/rownames you just identified above
  ddata$labels <- merge(ddata$labels, meta.x, by.x="label", by.y=matchcol)
  # plot it like you mean it
  ggplot(segment(ddata)) +
    geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +
    theme_dendro() +
    geom_text(data=label(ddata), aes_string(x='x', y='y', label=labels.colname, color=colors.colname, hjust=-0.1), size=6)+
    scale_color_brewer(type = "seq", palette = "Set1")+
    coord_flip() + scale_y_reverse(expand=c(0.2, 50)) +
    theme(axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank()) +
    ggtitle(title)
}