
scan1snps_HS <- function(chr, genoprobs, map, pheno, kinship, addcovar, query_func, start, end, intcovar = NULL){

    ### SNPs Association
    if(!is.null(intcovar)){
      
       pca.snp <- scan1snps(genoprobs = genoprobs, 
                            map = map,
                            pheno = pheno,
                            kinship = K[[chr]],
                            addcovar = covar,
                            query_func=query_variants,
                            chr = chr,
                            start= start,
                            end = end,
                            keep_all_snps=TRUE)
      
       pca.snp.int <- scan1snps(genoprobs = genoprobs, 
                                map = map,
                                pheno = pheno,
                                kinship = K[[chr]],
                                addcovar = covar,
                                intcovar = int_covar, 
                                query_func=query_variants,
                                chr = chr,
                                start = start,
                                end = end,
                                keep_all_snps=TRUE)
      
      
      
      pca.snp$lod <- pca.snp.int$lod - pca.snp$lod
      
      top <- top_snps(pca.snp$lod, pca.snp$snpinfo)
      top <- top[order(top$lod, decreasing = TRUE),]
      
      return(list('pca.snp' = pca.snp, 'top.snps' = top))
      
    }else{
      
      pca.snp <- scan1snps(genoprobs = genoprobs, 
                           map = map, 
                           pheno = pheno, 
                           kinship = K[[chr]],
                           addcovar = covar, 
                           query_func=query_variants,
                           chr = chr, 
                           start= start, 
                           end = end, 
                           keep_all_snps=TRUE)
      
      top <- top_snps(pca.snp$lod, pca.snp$snpinfo)
      top <- top[order(top$lod, decreasing = TRUE),]
      
      return(list('pca.snp' = pca.snp, 'top.snps' = top))
      
    }   
}
