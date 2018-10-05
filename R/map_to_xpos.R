# convert map to x-axis positions for plot_scan1
map_to_xpos <-
  function(map, gap)
  {
    if(length(map)==1) return(map[[1]])

    chr_range <- vapply(map, range, c(0,1), na.rm=TRUE)

    result <- map[[1]]-chr_range[1,1] + gap/2
    for(i in 2:length(map)) {
      result <- c(result,
                  map[[i]] - chr_range[1,i] + gap + max(result, na.rm=TRUE))
    }
    result
  }
