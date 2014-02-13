create_report<-function(path,samples,out,title,condition){
  require(Nozzle.R1)
  report<-""
  report <- newCustomReport( title );
  
  #report <- addTo( report, addTo( newSection( "My Introduction" ),
  #   newParagraph( "Hello World! This is a paragraph of text!" ) ) );
  #report <- addTo( report, addTo( newSection( "My Methods" ),
  #    method ) );
  FRfqc<- i_create_fastqc_figure(path,samples,out)
  Trqc<- i_create_mapping_stats(path,samples,out)
  FRfrccov<- i_create_gene_coverage(path,samples,out)
  Tcounts<- i_create_count_top(path,samples,out)
  FRcounts <- i_create_distribution_counts(path,samples,out,condition)
  
  
  report <- addTo( 
    report, addTo( newSection( "Quality of samples", class="results" ),
                   addTo( newSubSection( "Quality of base calling" ), FRfqc),
                   addTo( newSubSection( "Mapping stats" ), Trqc),
                   addTo( newSubSection("Gene coverage"), FRfrccov),    
                   addTo( newSubSection( "Top expressed genes" ), Tcounts),
                   addTo( newSubSection( "Counts distribution" ), FRcounts[[1]]),
                   addTo( newSubSection( "Normalized counts distribution" ), FRcounts[[2]]),
                   addTo( newSubSection("Cumulative expression"), FRcounts[[3]]),
                   addTo( newSubSection( "Samples similarity" ), FRcounts[[4]])
    )
    
    
  );
  
  return(report)
  
}

write_report<-function(report,outfile){
  writeReport( report, filename=outfile);
  
}


add_DE<-function(report,defile,condition){
  report <- addTo( 
    report,
    addTo(newSection("DE",class="results"),
          addTo(newSubSection("DE genes"),Tedger),
          addTo(newSubSection("Quality of DE analysis"),FRpval)
    )
  )
  
  return(report)
}