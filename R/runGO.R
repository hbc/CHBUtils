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

#' Clean and print results from enrichGO
#'
#' @param ego result object from enrichGO
#' @param limit integer limiting the number
#' of genes in each category. Bigger names
#' will give broad categories.
print_enrichGO = function(ego, limit = 100){
    
    require(knitr)
    if (.isvalid(ego)){
        idx = reduce_cp(ego$geneID, limit)
        return(kable(ego[idx, 1:7]))
    }
}

#' Run GO using clusterprofiler
#' 
#' @param genes character vector with list of genes
#' @param org OrgDb object (org.Mm.eg.db)
#' @param from character indicanting the name gene type.
#' Options at: (keytypes(org.Mm.eg.db))
#' @param ont ontology to do, BP, MF, CC
#' @param universe character vector with list of genes to use as universe
runGO = function(genes, org, from="SYMBOL", ont="BP", universe=NULL){
    require(clusterProfiler)
    require(AnnotationDbi)
    require(Rgraphviz)
    stopifnot(class(universe)=="character")
    stopifnot(class(genes)=="character")
    .accnum = convertIDs(genes, from, "ENTREZID", org, "useFirst")
    .accnumUNI = convertIDs(universe, from, "ENTREZID", org, "useFirst")
    stopifnot(length(universe)>length(genes))
    stopifnot(sum(!is.na(.accnum))>0)
    stopifnot(sum(!is.na(.accnumUNI))>0)
    ego <- enrichGO(gene = .accnum[!is.na(.accnum)], universe = .accnumUNI,
                    OrgDb = org, ont = ont, pAdjustMethod = "BH",
                    pvalueCutoff = 0.01, qvalueCutoff = 0.05, readable = TRUE)
    if("result" %in% slotNames(ego)){
        ego = simplify(ego)
        return(list(table=print_enrichGO(ego@result),
                    plot=dotplot(ego),
                    obj=ego))
        # plotGOgraph(ego)
    }
    
}
