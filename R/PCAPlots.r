##comment
PCAplot.eset <- function(eset=NULL, categories=NULL, title=NULL, colorpalette=NULL, alpha=1, numcomponents=4){
  alpha <- sprintf("%x", ceiling(alpha*255))
  colorpalette <- paste(colorpalette, alpha, sep="")
  eset.core <- exprs(eset) 
  myPca.core <- prcomp(t(eset.core))
  tmpPCAData.core <- as.data.frame(myPca.core$x[,1:numcomponents])
  pd <- pData(eset)
  colors <- colorpalette[factor(as.character(unlist(pd[,categories])))]
  legend_values=unique(cbind(colors, as.character(pd[,categories])))
  pairs(tmpPCAData.core, bg=colors, col="#606060", cex=2, pch=21, main=title, oma=c(8,5,5,14))
  legend("right", cex=0.7, col="#606060", pt.bg=legend_values[,1], pt.cex=1.5, legend=legend_values[,2],  pch=21, bty="n", x.intersp=1)
}
####################################
PCAplot.df <- function(df=NULL, meta.df=NULL, categories=NULL, title=NULL, colorpalette=NULL, alpha=1, numcomponents=6){
  alpha <- sprintf("%x", ceiling(alpha*255))
  colorpalette <- paste(colorpalette, alpha, sep="")
  myPca.core <- prcomp(t(df))
  tmpPCAData.core <- as.data.frame(myPca.core$x[,1:numcomponents])
  # SD of components
  colors <- colorpalette[factor(as.character(unlist(meta.df[,categories])))]
  legend_values=unique(cbind(colors, as.character(meta.df[,categories])))
  pairs(tmpPCAData.core, bg=colors, col="#606060", cex=2, pch=21, main=title, oma=c(8,5,5,14))
  legend("right", cex=0.7, col="#606060", pt.bg=legend_values[,1], pt.cex=1.5, legend=legend_values[,2],  pch=21, bty="n", x.intersp=1)
}

####################################
PCAplot.sd.eset <- function(eset=NULL,  title=NULL){
  eset.core <- exprs(eset)
  myPca.core <- prcomp(t(eset.core))
  # SD of components
  sdevdf <- data.frame(cbind(as.numeric(myPca.core$sdev),c(1:length(myPca.core$sdev))))
  sdevdf$prop <-  sdevdf$X1/sum(sdevdf$X1)
  sdevdf$cum <- cumsum(sdevdf$prop)
  ggplot(sdevdf, aes(x=X2, y=prop)) + 
    geom_point(size=4, color="red") + 
    scale_x_continuous('Component') + 
    scale_y_continuous('Standard Deviation') +
    ggtitle(title) +
    geom_line(data=sdevdf, aes(x=X2, y=cum))
}

####################################
PCAplot.sd.df <- function(eset=NULL, title=NULL){
  myPca.core <- prcomp(t(df))
  # SD of components
  sdevdf <- data.frame(cbind(as.numeric(myPca.core$sdev),c(1:length(myPca.core$sdev))))
  sdevdf$prop <-  sdevdf$X1/sum(sdevdf$X1)
  sdevdf$cum <- cumsum(sdevdf$prop)
  ggplot(sdevdf, aes(x=X2, y=prop)) + 
    geom_point(size=4, color="red") + 
    scale_x_continuous('Component') + 
    scale_y_continuous('Standard Deviation') +
    ggtitle(title) +
    geom_line(data=sdevdf, aes(x=X2, y=cum))
}
