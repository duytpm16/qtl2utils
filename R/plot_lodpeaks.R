plot_lodpeaks <- function(lod.peaks, map, window, slide, title = NULL){



        ### Counting number of LOD above threshold within a given window across each chromosome
        lod_df <- list()
        for(i in unique(names(map))){

            # Finding floor of minimum marker position and ceiling of maximum marker position
            min <- round_any(min(map[[i]]), 1, f = floor)
            max <- round_any(max(map[[i]]), 4, f = ceiling)

            # Creating x-axis scale. min to max with slide (or 'by' in seq function)
            x_axis_scale <- seq(min, max, slide)
            chr <- rep(i, length(x_axis_scale))

            # Getting LOD peaks from chromosome i
            sub <- subset(lod.peaks, qtl.chr == i)

            # Creating dataframe of counts of lod peaks above threshold
            count <- vector()
            pos <- ((x_axis_scale+window)+x_axis_scale)/2
            for(j in 1:length(pos)){
              count[j] <- sum((abs(sub$qtl.pos - pos[j])  <= window/2))
            }

            lod_df[[i]] <- data.frame(chr = chr, pos = pos, count = count)
        }


        lod_df <- rbindlist(lod_df)
        lod_df$chr <- factor(lod_df$chr, levels = c(1:19,'X'))




        return(lod_df)

}
