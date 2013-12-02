####################################
col2rownames <- function(df, colname, removecol=FALSE){
  row.names(df) <- df[,colname]
  if(removecol){df[,colname] <- NULL}
  return(df)
}
