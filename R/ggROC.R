ggROC = function(df, facetName, groupName = "grp", predName = "res") {
  require(plyr)
  require(ggplot2)
  df = df[complete.cases(df),]
  plotdata = ddply(df, facetName,
      function(x) data.frame(roc(x[,groupName], x[,predName])[c("sensitivities",
                                                            "specificities")]))
  plotdata$specificities = 1 - plotdata$specificities
  colnames(plotdata) = c(facetName, "tpr", "fpr")
  p = ggplot(plotdata, aes(fpr, tpr)) + 
      geom_line(aes_string(colour=facetName))
  return(p)
}

