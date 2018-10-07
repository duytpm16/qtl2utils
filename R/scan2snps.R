
scan2snps <- function(genoprobs, map, pheno, kinship = NULL, addcovar = NULL,
   Xcovar = NULL, intcovar = NULL, weights = NULL, reml = TRUE,
   model = c("normal", "binary"), query_func = NULL, chr = NULL,
   start = NULL, end = NULL, snpinfo = NULL, batch_length = 20,
   keep_all_snps = FALSE, cores = 1, ...){


    ### SNPs Association
    if(!is.null(intcovar)){

       pca.snp <- scan1snps(genoprobs, map, pheno, kinship = NULL, addcovar = NULL,
                            Xcovar = NULL, intcovar = NULL, weights = NULL, reml = TRUE,
                            model = c("normal", "binary"), query_func = NULL, chr = NULL,
                            start = NULL, end = NULL, snpinfo = NULL, batch_length = 20,
                            keep_all_snps = FALSE, cores = 1, ...)

       pca.snp.int <- scan1snps(genoprobs, map, pheno, kinship = NULL, addcovar = NULL,
                                Xcovar = NULL, intcovar = NULL, weights = NULL, reml = TRUE,
                                model = c("normal", "binary"), query_func = NULL, chr = NULL,
                                start = NULL, end = NULL, snpinfo = NULL, batch_length = 20,
                                keep_all_snps = FALSE, cores = 1, ...)



      pca.snp$lod <- pca.snp.int$lod - pca.snp$lod

      top <- top_snps(pca.snp$lod, pca.snp$snpinfo)
      top <- top[order(top$lod, decreasing = TRUE),]

      return(list('pca.snp' = pca.snp, 'top.snps' = top))

    }else{

      pca.snp <- scan1snps(genoprobs, map, pheno, kinship = NULL, addcovar = NULL,
                           Xcovar = NULL, intcovar = NULL, weights = NULL, reml = TRUE,
                           model = c("normal", "binary"), query_func = NULL, chr = NULL,
                           start = NULL, end = NULL, snpinfo = NULL, batch_length = 20,
                           keep_all_snps = FALSE, cores = 1, ...)

      top <- top_snps(pca.snp$lod, pca.snp$snpinfo)
      top <- top[order(top$lod, decreasing = TRUE),]

      return(list('pca.snp' = pca.snp, 'top.snps' = top))

    }
}
