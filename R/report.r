i_create_fastqc_figure<-function(path,samples,out){
  require(ggplot2)
  require(Nozzle.R1)
  ############QUALITY PER NT###############
  tfqc<-data.frame()
  for (f in samples){
    src<-paste(path,f,"/qc/fastqc/fastqc_data.txt",sep="")
    endnt<-system(paste("grep -n END ",src," | head -2 | tail -1| sed 's/:/\t/' | cut -f 1"),
                  intern = TRUE)
    endnt<-as.integer(endnt)
    fqc<-suppressMessages(read.table(src,nrows=endnt-14,skip=12,sep="\t"))
    names(fqc)<-c("nt","mean","median","lowQ","upQ","10Q","90Q")
    fqc$sample<-substr(f,1,21)
    tfqc<-rbind(tfqc,fqc)
    
  }
  tfqcFile<-paste(out,"tfqc.txt",sep="")
  write.table(tfqc,tfqcFile,row.names=F,quote=F,sep="\t")
  
  
  p<-ggplot(tfqc,aes(factor(nt),mean)) +
    #geom_boxplot() + 
    geom_jitter(position=position_jitter(width=0.3),aes(factor(nt), mean,colour=factor(sample)))+
    geom_abline(intercept = 30,slope=0)+
    ylim(25,42)+
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_blank())+
    labs(list(x="nucleotide",y="score",colour="Samples"))
  
  
  ffqcFile="ffqc.jpg"
  ffqcHFile="ffqcH.pdf"
  jpeg(paste(out, ffqcFile, sep="" ) ,width=600,height=400,quality=100 );
  print(p)
  dev.off();
  pdf(paste( out, ffqcHFile, sep="" )  );
  print(p)
  dev.off();
  
  # create a figure and make it available for exporting
  FRfqc <- newFigure( ffqcFile, fileHighRes=ffqcHFile, exportId="FIGURE_FASTQC",
                      "This figure shows the quality per nucleotide for each sample.
                      Any value above 30 is quite good quality." );
  ################################################
  return(FRfqc)
}

i_create_rnaseqc_figure<-function(path,samples,out){
  require(ggplot2)
  require(Nozzle.R1)
  
  ############COVERAGE###############
  
  trqccov<-data.frame()
  for (f in samples){
    print (f)
    src<-paste(path,f,"/qc/rnaseqc/meanCoverageNorm_low.txt",sep="")
    temp<-read.table(src,skip=1)
    temp$type<-"low"
    temp$pos<-1:nrow(temp)
    temp$sample<-substr(f,1,21)
    trqccov<-rbind(trqccov,temp)
    src<-paste(path,f,"/qc/rnaseqc/meanCoverageNorm_medium.txt",sep="")
    temp<-read.table(src,skip=1)
    temp$type<-"medium"
    temp$pos<-1:nrow(temp)
    temp$sample<-substr(f,1,21)
    trqccov<-rbind(trqccov,temp)
    src<-paste(path,f,"/qc/rnaseqc/meanCoverageNorm_high.txt",sep="")
    temp<-read.table(src,skip=1)
    temp$type<-"high"
    temp$pos<-1:nrow(temp)
    temp$sample<-substr(f,1,21)
    trqccov<-rbind(trqccov,temp)
    
  }
  
  trqccovFile<-paste(out,"trqccov.txt",sep="")
  write.table(trqccov,trqccovFile,row.names=F,quote=F,sep="\t")
  
  trqccov$type<-factor(trqccov$type,levels=c("low","medium","high"))
  p<-ggplot(trqccov,aes(pos,log2(V1+0.5),colour=sample)) +
    geom_line()+
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_blank())+
    labs(list(x="gene",y="coverage",colour="Samples"))+
    facet_wrap(~type,nrow=3)
  
  
  stdFile="frqccov.jpg"
  highFile= "frqccovH.jpg"
  jpeg( paste( out,stdFile , sep="" ),width=600,height=400,quality=100 );
  print(p)
  dev.off(frqccovFile);
  pdf( paste( out, highFile, sep="" ));
  print(p)
  dev.off();
  
  # create a figure and make it available for exporting
  FIG <- newFigure( stdFile, fileHighRes=highFile,exportId="FIGURE_COV",
                         "This shows the gene coverage by reads. 
                    A good signal is an equal coverage along gene." );
  ################################################
  return(FIG)
  
}

i_create_gene_coverage<-function(path,samples,out){
  require(ggplot2)
  require(Nozzle.R1)
  
  ############COVERAGE###############
  trqccov<-data.frame()
  for (f in samples){
    src<-paste(path,f,"/qc/rnaseqc/meanCoverageNorm_low.txt",sep="")
    temp<-read.table(src,skip=1)
    temp$type<-"low"
    temp$pos<-1:nrow(temp)
    temp$sample<-substr(f,1,21)
    trqccov<-rbind(trqccov,temp)
    src<-paste(path,f,"/qc/rnaseqc/meanCoverageNorm_medium.txt",sep="")
    temp<-read.table(src,skip=1)
    temp$type<-"medium"
    temp$pos<-1:nrow(temp)
    temp$sample<-substr(f,1,21)
    trqccov<-rbind(trqccov,temp)
    src<-paste(path,f,"/qc/rnaseqc/meanCoverageNorm_high.txt",sep="")
    temp<-read.table(src,skip=1)
    temp$type<-"high"
    temp$pos<-1:nrow(temp)
    temp$sample<-substr(f,1,21)
    trqccov<-rbind(trqccov,temp)
    
  }
  
  trqccovFile<-paste(out,"trqccov.txt",sep="")
  write.table(trqccov,trqccovFile,row.names=F,quote=F,sep="\t")
  
  trqccov$type<-factor(trqccov$type,levels=c("low","medium","high"))
  p<-ggplot(trqccov,aes(pos,log2(V1+0.5),colour=sample)) +
    geom_line()+
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_blank())+
    labs(list(x="gene",y="coverage",colour="Samples"))+
    facet_wrap(~type,nrow=3)
  
  
  stdFile="frqccov.jpg"
  highFile="frqccovH.pdf"
  jpeg( paste( out ,stdFile , sep="" ) ,width=600,height=400,quality=100 );
  print(p)
  dev.off();
  pdf( paste( out , highFile, sep="" ) );
  print(p)
  dev.off();
  
  # create a figure and make it available for exporting
  FIG <- newFigure( stdFile, fileHighRes=highFile,exportId="FIGURE_COV",
                         "This shows the gene coverage by reads. 
                    A good signal is an equal coverage along the gene." );
  ################################################
  
  return(FIG)
}

i_create_count_top<-function(path,ssamples,out){
  require(Nozzle.R1)
  
  mainfolder<-dir(path,pattern="project")
  ############COUNTS###############
  counts<-read.table(paste(sep="",path,mainfolder,"/annotated_combined.counts"),header=T,sep="\t")
  names(counts)<-substr(names(counts),1,21)
  counts$mean<-rowMeans(counts[,2:(ncol(counts)-1)])
  counts<-counts[order(-counts$mean),]
  counts<-counts[counts$mean>0,]
  
  countsFile<-"counts.txt"
  write.table(counts,paste(out,countsFile,sep=""),row.names=F,quote=F,sep="\t")
  short<-counts[1:30,]
  TAB <- newTable(short , file=countsFile, exportId="TABLE_COUNTS",
                      "Top genes" );
  
  return(TAB)
}

i_create_distribution_counts<-function(path,samples,out,condition){
  require(ggplot2)
  require(Nozzle.R1)
  require(reshape)
  require(DESeq2)
  
  mainfolder<-dir(path,pattern="project")
  
  counts<-read.table(paste(sep="",path,mainfolder,"/combined.counts"),header=T,row.names=1,sep="\t")
  names(counts)<-substr(names(counts),1,21)
  keep<-rowSums(counts>0)>=1
  counts<-counts[keep,]
  cshape<-suppressMessages(melt(counts,by=1))
  p<-ggplot(cshape,aes(factor(variable),log2(value+0.1))) +
    geom_boxplot() + 
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_blank()) +
    labs(list(x="",y="log2(raw counts)"))
  
  stdFile= "cbox.jpg"
  jpeg( paste( out,stdFile, sep="" ),width=600,height=400,quality=100 );
  print(p)
  dev.off();
  
  # create a figure and make it available for exporting
  FIG1 <- newFigure( stdFile,exportId="FIGURE_COUNTS",
                    "Distribution of raw counts" );
  
  ma<-counts
  design<-data.frame(row.names=names(ma),condition=condition)
  dse<-DESeqDataSetFromMatrix(countData = ma,
                              design=~condition,
                              colData=design)
  rld<-rlogTransformation(dse)
  p<-ggplot(suppressMessages(melt(assay(rld))),aes(factor(X2),value)) +
    geom_boxplot() + 
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_blank())+
    labs(list(x="",y="log2(normalized counts)"))
  
  stdFile= "cboxnorm.jpg"
  jpeg( paste( out,stdFile, sep="" ),width=600,height=400,quality=100 );
  print(p)
  dev.off();
  
  # create a figure and make it available for exporting
  FIG2 <- newFigure( stdFile,exportId="FIGURE_COUNTSNORM",
                           "Distribution of normalized counts" );
  
  cum<-assay(rld)
  cum<-as.data.frame(apply(cum,2,sort))
  cum$pos<-nrow(cum):1
  cum<- suppressMessages(melt(cum,id.vars="pos"))
  p<-ggplot(cum,aes(pos,value,colour=variable)) +
    geom_point() +
    theme_bw(base_size = 12) +
    labs(list(x="sorted genes",y="log2(normalized counts)"))
  
  sdtFile="fcum.jpg"
  highFile="fcumH.pdf"
  jpeg(paste( out,sdtFile , sep="" ) ,width=600,height=400,quality=100 );
  print(p)
  dev.off()
  pdf( paste( out, highFile, sep="" ) )
  print(p)
  dev.off()
  
  # create a figure and make it available for exporting
  FIG3<- newFigure( sdtFile, fileHighRes=highFile,exportId="FIGURE_CUM",
                       "Cumulative expression by sample." )
  

  p<-mds(assay(rld),condition)
  p<-p + theme_bw(base_size = 12) +
    geom_text(aes(one,two,label=label))+
    scale_color_brewer(palette="Set1")
  
  stdFile= "mds.jpg"
  jpeg( paste( out,stdFile, sep="" ),width=600,height=400,quality=100 )
  print(p)
  dev.off()
  
  # create a figure and make it available for exporting
  FIG4 <- newFigure( stdFile,exportId="FIGURE_MDS",
                       "Multidimensional scaling plot" );
  
  
  return(list(FIG1,FIG2,FIG3,FIG4))
  
}

i_create_mapping_stats<-function(path,samples,out){
  require(Nozzle.R1)
  
  ############ALIGNED#################
  trqc<-data.frame()
  for (f in samples){
    src<-paste(path,f,"/qc/rnaseqc/metrics.tsv",sep="")
    rqc<-read.table(src,sep="\t",skip=1)
    rqc$sample<-substr(f,1,21)
    trqc<-rbind(trqc,rqc)
  }
  
  trqc<-trqc[,c(1,19,17,8,20,23,7,38)]
  names(trqc)<-c("Sample","Total","Mapped","Mapping rate","rRNA","transcripts","Exonic rate","Bases MM rate")
  
  tFile<-"trqc.txt"
  write.table(trqc,paste(out,tFile,sep=""),row.names=F,quote=F,sep="\t")
  
  # create a table and make it available for exporting
  TAB <- newTable( trqc, file=tFile, exportId="TABLE_RNASEQC",
                    "Stats of mapping reads " );
  
  
  ################################################
  
  return(TAB)
  
}


