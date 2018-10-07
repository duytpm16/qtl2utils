

scan2blup <- function(chr, marker.id, genoprobs, pheno, kinship = NULL, addcovar = NULL,
                      intcovar = NULL, samples = NULL, covar.factors = NULL, nullcovar = NULL, contrasts = NULL, se = FALSE, reml = TRUE,
                      tol = 0.000000000001, cores = 1, quiet = TRUE){
### Alle Effects
if(!is.null(intcovar)){

    stopifnot(!is.null(samples) & !is.null(covar.factors) & !is.null(covar_name))

    male_mouse <- samples[samples$sex == 'M',]
    female_mouse <- samples[samples$sex == 'F',]

    f <- as.formula(paste0('~',paste0(covar.factors[!(covar.factors %in% int_covar_name)], collapse = '+')))

    male_addcovar <- model.matrix(f, data = male_mouse)[,-1]
    female_addcovar <- model.matrix(f, data = female_mouse)[,-1]

    pheno_male <- pheno[rownames(pheno) %in% male_mouse$mouse.id,, drop = FALSE]
    pheno_female <- pheno[rownames(pheno) %in% female_mouse$mouse.id,, drop = FALSE]


    m_blup = scan1blup(genoprobs = genoprobs[,chr], pheno = pheno_male,
                       kinship = K[[chr]], addcovar = male_addcovar, cores = cores)

    f_blup = scan1blup(genoprobs = genoprobs[,chr], pheno = pheno_female,
                       kinship = K[[chr]], addcovar = female_addcovar, cores = cores)


    gp = genoprobs[,chr]
    gp[[1]] = gp[[1]][,,marker.id,drop=FALSE]
    coeff_blup_m = scan1blup(genoprobs = gp,
                             pheno = pheno_male,
                             kinship = K[[chr]],
                             addcovar = male_addcovar,
                             cores = num_cores)


    coeff_blup_f = scan1blup(genoprobs = gp,
                             pheno = pheno_female,
                             kinship = K[[chr]],
                             addcovar = female_addcovar,
                             cores = num_cores)

    return(list('male_blupscan' = m_blup, 'male_blup_coeff' = coeff_blup_m, 'female_blupscan' = f_blup, 'female_blup_coeff' = coeff_blup_f))

}else{

    # Do allele effects for additive scan only
    pca_blup <- scan1blup(genoprobs = genoprobs[,chr],
                          pheno = pheno,
                          kinship = K[[chr]],
                          addcovar = covar,
                          cores = num_cores)


    gp = genoprobs[,chr]
    gp[[1]] = gp[[1]][,,marker.id,drop=FALSE]
    blup = scan1blup(genoprobs = gp,
                     pheno = pheno,
                     kinship = K[[chr]],
                     addcovar = covar,
                     cores = num_cores)

    return(list('blup_scan' = pca_blup, 'blup_coeff' = blup))

}
}
