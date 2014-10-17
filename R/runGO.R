#' Prepare universe for runGO
#'
#' @param anno_db GO annotation of the organism: org.Dr.egGO
#' @param organism Name of the organism: Denio reria
#' @export
#' #' @examples
#' \dontrun{
#' library(org.Dr.eg.db)
#' g<-read_entrez_id()
#' set <- universeGO(org.Dr.egGO,"Denio reiro")
#' }

universeGO<-function(anno_db,organims)
{
  suppressMessages(require(Category))
  suppressMessages(require(GOstats))
  suppressMessages(require(GSEABase))
  
  frame = toTable(anno_db)
  goframeData = data.frame(frame$go_id, frame$Evidence, frame$gene_id)
  goFrame = GOFrame(goframeData,organism=organims)
  goAllFrame = GOAllFrame(goFrame)
  
  gsc <- GeneSetCollection(goAllFrame, setType = GOCollection())
  universe <- Lkeys(anno_db)  
  list(gsc,universe)
}

#' Run basic GO enrichment of model organism
#'
#' @param g list of entrezid genes 
#' @param set list from create_universe
#' @param gotree list from as.list(GOMFCHILDREN). This will determine if
#' the more generic or specific term is kept.
#' @param type MF, BP or CC
#' @return data.frame with significant GO terms 
#' @export
#' @examples
#' \dontrun{
#' rwa <- runGO(g.ezid,set,"MF",as.list(GOMFCHILDREN))
#' }
runGO<-function(g,set,type,gotree)
{
  gsc <- set[[1]]
  universe <- set[[2]]
  go <- do_go(g,gsc,universe,type)
  list(reduce(g,go,gotree),go[[2]])
}

do_go <- function(g,gsc,universe,type)
{
  params <- suppressWarnings(GSEAGOHyperGParams(name="go",
                                                geneSetCollection=gsc,
                                                geneIds = as.character(g),
                                                universeGeneIds = universe,
                                                ontology = type,
                                                pvalueCutoff = 0.05,
                                                conditional = FALSE,
                                                testDirection = "over"))
  
  Over <- hyperGTest(params)
  tab <- GOstats::summary(Over)
  tab$fdr <- p.adjust(tab$Pvalue, method="BH")
  list(tab,Over)
}

reduce <- function(g,go,childs)
{
  res <- data.frame()
  f <- subset(go[[1]],fdr<0.1 & Count > 5 & OddsRatio > 1.5 )
  for (go_id in f[,1])
  {
    if ( sum(unlist(childs[go_id]) %in% f[,1])<=1){
      info <- f[f[,1]==go_id,]
      genes_in<-go[[2]]@goDag@nodeData@data[go_id][[1]][1][[1]]
      gene_in_cit <- intersect(as.character(g),genes_in)
      res <- rbind(res,data.frame(
        go = go_id,
        term = info$Term,
        OddsRatio = info$OddsRatio,
        FDR = info$fdr
      ))
    }
  }
  res
}


plotGO <- function(go, name, ext = 'png'){
    suppressMessages(require(RamiGO))
    t <- getAmigoTree (goIDs = go[[1]]$go, pvalues = go[[1]]$FDR,
                             pcolors = c('white', 'magenta'),
                             psplit = c(0.1, 0.05, 0.025, 0.01, 0.00025),
                             filename = name, picType = ext,
                             modeType = 'amigo', saveResult = TRUE)
}

#' Run basic GO enrichment of model organism
#'
#' @param go object from runGO
#' @param gene character with string to match
#' @return data.frame with matched terms
#' @export
#' @examples
#' \dontrun{
#' getGO(go,"hemato")
#' }
getGO <- function(go,gene){
    go[[1]][grepl(gene,go[[1]]$term),]
}