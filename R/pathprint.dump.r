library(pathprint)
library(reshape2)
library(biomaRt)
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")

paths.l <- pathprint.Hs.gs
temp <- llply(paths.l, function(n) {
  getBM(attributes="hgnc_symbol", filters="entrezgene", values=n, mart=ensembl)
})


rowMax <- max(sapply(temp, length)) 
temp2 <- do.call(rbind, lapply(temp, function(x){ 
  length(x) <- rowMax 
  x 
}))

write.table(temp2, file="pathprint.pathways.Hs.tab", sep="\t", quote=F, row.names=T, col.names=F)
