####################################
rownames2col <- function(df, colname) {
  output <- cbind(row.names(df), df)
  colnames(output)[1] <- colname
  return(output)
}
