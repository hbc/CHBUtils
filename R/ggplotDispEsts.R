#' Plot Dispersion Estimates - nicer ggplot version of DESeq graph
#'
#' @param cds CountDataSet, required
#' @return plot
#' @seealso \code{\link{DESeq}} for description of the Dispersion plots
#' @export
#' @examples
#' ggplotDispEsts(cds)


ggplotDispEsts = function(cds) {
  estimates = data.frame(means = rowMeans(counts(cds, normalized=TRUE)),
                         variance = fitInfo(cds)$perGeneDispEsts)
  xg = 10^seq(-0.5, 5, length.out=300)
  yg = fitInfo(cds)$dispFun(xg)
  fitline = data.frame(xg=xg, yg=yg)
  p = ggplot(estimates, aes(means, variance)) + geom_point(size=1, alpha=0.4) +
    scale_x_log10() + scale_y_log10() +
    geom_line(data=fitline, aes(xg, yg), color="red") +
    labs(title="Dispersion estimation while pooling all samples") +
    xlab("mean number of mapped reads per gene") +
    ylab("estimated dispersion")
  p
}