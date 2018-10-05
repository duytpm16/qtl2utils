
scan1_HS <- function(genoprobs, pheno, kinship, addcovar, intcovar = NULL, cores = 1){
  
        ### Genome scan
        if(!is.null(intcovar)){
           # Do difference between interaction and additive scan
           pca.scan <- scan1(genoprobs = genoprobs, 
                             pheno = pheno, 
                             kinship = K, 
                             addcovar = covar, 
                             cores = num_cores)
          
          pca.scan.int <- scan1(genoprobs = genoprobs, 
                                pheno = pheno, 
                                kinship = K, 
                                addcovar = covar,
                                intcovar = int_covar,
                                cores = num_cores)
          
          pca.scan.diff <- pca.scan.int - pca.scan
          
          return(list('add.scan' = pca.scan, 'int.scan' = pca.scan.int, 'diff.scan' = pca.scan.diff))
          
        }else{
          # Do additive scan only
          pca.scan <- scan1(genoprobs = genoprobs, pheno = pheno, kinship = K, addcovar = covar, cores = num_cores)
          
          return(pca.scan)
        }
}
