#' add a column to a dataframe consisting of the current rownames
#'
#' @param df dataframe
#' @param colname character vector specifying id of column with new rownames
#' @param removecol boolean specifying wheter to remove the column used to rename rownames
#' @return dataframe with column contents as rownames
#' @export
#' @examples
#' rownamess2col(dataframe, "newcolumnid")


####################################
col2rownames <- function(df, colname, removecol=TRUE) {
  output <- df
  row.names(output) <- output[,colname]
  if (removecol) {
    output[,colname] <- NULL
  }
  return(output)
}
