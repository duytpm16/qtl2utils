scan1_snpasso_rng <- function(markers, scan1_output, max_chr, start, end){
   
    new_markers <- subset(markers, chr == max_chr & (pos >= start & pos <= end))
    
    new_markers <- transform(merge(new_markers, scan1_output,by=0,all.x=TRUE, all.y = FALSE), 
                         row.names=pos, Row.names=NULL)
    
    scan1_rng <- new_markers[, !colnames(new_markers) %in% colnames(markers), drop = FALSE]
    scan1_rng$pos <- as.numeric(rownames(scan1_rng))

    
    
    return(scan1_rng)
  
}
