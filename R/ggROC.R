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

jaccard = function(df, column, cutoff=0.05) {
    algorithms = unique(df$algorithm)
    sig = list()
    for(alg in algorithms) {
        sig[alg] = list(subset(df, algorithm == alg & df[,column] < cutoff)$id)
    }
    ji = list()
    for(alg in algorithms) {
        for(alg2 in algorithms) {
            if (alg == alg2) next
            algs = sort(c(alg, alg2))
            alg1 = algs[1]
            alg2 = algs[2]
            comparison = paste(alg1, "_vs_", alg2, sep="")
            alg_intersection = intersect(sig[alg1][[1]], sig[alg2][[1]])
            alg_union = union(sig[alg1][[1]], sig[alg2][[1]])
            ji[comparison] = length(alg_intersection) / length(alg_union)
        }
    }
    return(data.frame(ji))
}

ji = jaccard(df, "padj")
