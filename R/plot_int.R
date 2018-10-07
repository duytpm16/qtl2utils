
plot_int <- function(add_scan, int_scan, diff_scan, map, lodcolumn = 1, chr = NULL){

      ### Plot interaction, additive, and difference scans
      plot_scan1(x = int_scan,
                 map = map,
                 lodcolumn = lodcolumn,
                 chr = chr,
                 main = paste('Chromsome',chr,'Scan'),
                 col = 'maroon',
                 ylim = c(0, max(int_scan[,lodcolumn]) + 5))

      plot_scan1(x = add_scan,
                 add = TRUE,
                 map = map,
                 lodcolumn = lodcolumn,
                 chr = chr,
                 main = paste('Chromsome',chr,'Scan'))



      plot_scan1(x = diff_scan,
                 add = TRUE,
                 map = map,
                 lodcolumn = lodcolumn,
                 chr = chr,
                 main = paste('Chromsome',chr,'Scan'),
                 col = '#66CC00')



      ### Plotting legend
      pos <- max_scan1(diff_scan, map = map, lodcolumn = lodcolumn, chr = chr, na.rm = TRUE)$pos
      chr_range <- range(map[[chr]])
      chr_seg <- c(chr_range[1], (chr_range[2]+chr_range[1])/2, chr_range[2], pos)
      seg_dist <- as.matrix(dist(chr_seg))[,-4]
      min_dist_to_pos <- which.min(seg_dist[4,])
      if(min_dist_to_pos == 2 | min_dist_to_pos == 1){
        legend("topright", inset = .05, legend = c('Sex Interaction','Difference','Additive'),
               col = c('maroon','#66CC00','darkslateblue'), lty = 1, lwd = 3,
               bg = 'white', box.col = 'white')
      }else{
        legend("topleft", inset = .05, legend = c('Sex Interaction','Difference','Additive'),
               col = c('maroon','#66CC00','darkslateblue'), lty = 1, lwd = 3,
               bg = 'white', box.col = 'white')
      }
}
