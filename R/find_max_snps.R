
find_max_snps <- function(snp_pos, snp_id, snp_lod, max_peak = NULL, haplo_peak_drop = NULL){
  
    snp_pos <- data.frame(pos = snp_pos)
    rownames(snp_pos) <- snp_id
    
    snp_data <- merge(snp_pos, snp_lod, by = 'row.names')
    colnames(snp_data)[ncol(snp_data)] <- 'snp_lod'
    
    if(!is.null(haplo_peak_drop)){
       max_snp <- snp_data[which.max(snp_data[,'snp_lod']),]
       
       
       max_snp_peak_drop <- subset(snp_pos, pos >= haplo_peak_drop[1] & pos <= haplo_peak_drop[2]) 
       max_snp_peak_drop <- max_snp_peak_drop[which.max(max_snp_peak_drop[,3]),]
       
       
       all_max_snp <- rbind(max_snp,max_snp_peak_drop)
       return(all_max_snp)
    }else{
       max_snp <- snp_data[which.max(snp_data[,'snp_lod']),]
       return(max_snp)
       
    }
  
}

