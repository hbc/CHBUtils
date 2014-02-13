#' add a column to a dataframe consisting of the current rownames
#'
#' @param df dataframe
#' @param colname character vector specifying new column id
#' @return dataframe with rownames as a column
#' @export
#' @examples
#' rownamess2col(dataframe, "newcolumnid")


####################################
row2colnames <- function(df, colname) {
  output <- cbind(row.names(df), df)
  colnames(output)[1] <- colname
  return(output)
}
