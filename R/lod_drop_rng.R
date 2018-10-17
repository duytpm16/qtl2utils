lod_drop_rng <- function(markers, scan1_output, max_chr, max_peak, hilit_drop){

    markers <- transform(merge(markers[,c('chr','pos')], scan1_output ,by=0,all=TRUE),
                          row.names=Row.names, Row.names=NULL)

    markers <- subset(markers, chr == max_chr)
    lod_drop <- max_peak - 1.5

    markers <- markers[markers[,ncol(markers)] >= lod_drop,]

    return(range(markers$pos))

}
