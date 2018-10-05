plot_scan2 <-
  function(x, map, snp_lod, snpinfo, max_peak, haplo_peak_drop = NULL, lodcolumn=1, chr=NULL, add=FALSE, gap=25, ...)
  {
    if(is.null(map)) stop("map is NULL")

    if(!is.list(map)) map <- list(" "=map) # if a vector, treat it as a list with no names

    # subset chromosomes
    if(!is.null(chr)) {
      chri <- match(chr, names(map))
      if(any(is.na(chri)))
        stop("Chromosomes ", paste(chr[is.na(chri)], collapse=", "), " not found")
      x <- subset_scan1(x, map, chr)
      map <- map[chri]
    }


    # align scan1 output and map
    tmp <- align_scan1_map(x, map)
    x <- tmp$scan1
    map <- tmp$map

    if(nrow(x) != length(unlist(map)))
      stop("nrow(x) [", nrow(x), "] != number of positions in map [",
           length(unlist(map)), "]")

    # pull out lod scores
    if(length(lodcolumn)==0) stop("lodcolumn has length 0")
    if(length(lodcolumn) > 1) { # If length > 1, take first value
      warning("lodcolumn should have length 1; only first element used.")
      lodcolumn <- lodcolumn[1]
    }
    if(is.character(lodcolumn)) { # turn column name into integer
      tmp <- match(lodcolumn, colnames(x))
      if(is.na(tmp)) stop('lodcolumn "', lodcolumn, '" not found')
      lodcolumn <- tmp
    }
    if(lodcolumn < 1 || lodcolumn > ncol(x))
      stop("lodcolumn [", lodcolumn, "] out of range (should be in 1, ..., ", ncol(x), ")")
    lod <- unclass(x)[,lodcolumn]

    # internal function; trick to be able to pull things out of "..."
    #    but still have some defaults for them
    plot_scan2_internal <-
      function(snp_lod, snpinfo, map, lod, max_peak, haplo_peak_drop = NULL,
               add=FALSE, gap,bgcolor="gray90", altbgcolor="gray85",lwd=2,
               col="darkslateblue", altcol=NULL, xlab=NULL, ylab="LOD score",
               xlim=NULL, ylim=NULL, xaxs="i", yaxs="i", main="",
               mgp.x=c(2.6, 0.5, 0), mgp.y=c(2.6, 0.5, 0), mgp=NULL, las=1,
               hlines=NULL, hlines_col="white", hlines_lwd=1, hlines_lty=1,
               vlines=NULL, vlines_col="white", vlines_lwd=1, vlines_lty=1,
               sub="", ...)
      {
        dots <- list(...)
        onechr <- (length(map)==1) # single chromosome

        xpos <- map_to_xpos(map, gap)
        chrbound <- map_to_boundaries(map, gap)

        if(!add) { # new plot
          if(is.null(ylim))
            ylim <- c(0, max(lod, na.rm=TRUE)*1.02)

          if(is.null(xlim)) {
            xlim <- range(xpos, na.rm=TRUE)
            if(!onechr) xlim <- xlim + c(-gap/2, gap/2)
          }

          if(is.null(xlab)) {
            if(onechr) {
              if(names(map) == " ") xlab <- "Position"
              else xlab <- paste("Chr", names(map), "position")
            }
            else xlab <- "Chromosome"
          }

          # margin parameters
          if(!is.null(mgp)) mgp.x <- mgp.y <- mgp
          #View(data.frame(xpos = xpos, lod = lod))
          # make basic plot

          plot(xpos, lod, xlab="", ylab="", xlim=xlim, ylim=ylim,
               xaxs=xaxs, yaxs=yaxs, xaxt="n", yaxt="n", type="p",
               main=main, sub=sub)


          max_snps <- find_max_snps(snp_pos = snpinfo$pos, snp_id = snpinfo$snp_id, snp_lod = snp_lod, max_peak = max_peak, haplo_peak_drop = NULL)

          snp_x <- c(max_snps$pos)
          snp_y <- c(max_snps$snp_lod)


          # add background rectangles
          u <- par("usr")
          if(!is.null(bgcolor))
            rect(u[1], u[3], u[2], u[4], col=bgcolor, border=NA)
          if(!is.null(altbgcolor) && !onechr) {
            for(i in seq(2, ncol(chrbound), by=2))
              rect(chrbound[1,i], u[3], chrbound[2,i], u[4], col=altbgcolor, border=NA)
          }

          # include axis labels?
          if(is.null(dots$xaxt)) dots$xaxt <- par("xaxt")
          if(is.null(dots$yaxt)) dots$yaxt <- par("yaxt")

          # add x axis unless par(xaxt="n")
          if(dots$xaxt != "n") {
            if(onechr) {
              axis(side=1, at=pretty(xlim), mgp=mgp.x, las=las, tick=FALSE)
            }
            else {
              loc <- colMeans(chrbound)
              odd <- seq(1, length(map), by=2)
              even <- seq(2, length(map), by=2)
              axis(side=1, at=loc[odd], names(map)[odd],
                   mgp=mgp.x, las=las, tick=FALSE)
              axis(side=1, at=loc[even], names(map)[even],
                   mgp=mgp.x, las=las, tick=FALSE)
            }
          }

          # add y axis unless par(yaxt="n")
          if(dots$yaxt != "n") {
            axis(side=2, at=pretty(ylim), mgp=mgp.y, las=las, tick=FALSE)
          }
          dots$xaxt <- dots$yaxt <- NULL # delete those

          # x and y axis labels
          title(xlab=xlab, mgp=mgp.x)
          title(ylab=ylab, mgp=mgp.y)


          # grid lines
          if(onechr && !(length(vlines)==1 && is.na(vlines))) { # if vlines==NA (or mult chr), skip lines
            if(is.null(vlines)) vlines <- pretty(xlim)
            abline(v=vlines, col=vlines_col, lwd=vlines_lwd, lty=vlines_lty)
            abline(v = max_peak)
            if(!is.null(haplo_peak_drop)){
              abline(v = haplo_peak_drop[1], lty = 3)
              abline(v = haplo_peak_drop[2], lty = 3)
            }
            segments(x0 = snp_x, y0= 0,x1 = snp_x, y1 =snp_y, col = 'violetred')
            text_pos <- as.matrix(dist(c(haplo_peak_drop[1],haplo_peak_drop[2], max_peak)))[,-3]

          }
          if(!(length(hlines)==1 && is.na(hlines))) { # if hlines==NA, skip lines
            if(is.null(hlines)) hlines <- pretty(ylim)
            abline(h=hlines, col=hlines_col, lwd=hlines_lwd, lty=hlines_lty)
            abline(v = max_peak)
            segments(x0 = snp_x, y0= rep(0,length(snp_x)),x1 = snp_x, y1 =snp_y, col = 'violetred')

          }

        }

        # plot each chromosome
        indexes <- map_to_index(map)
        for(i in seq(along=indexes)) {
          # if altcol provided, have chromosomes alternate colors
          if(is.null(altcol) || i %% 2) this_col <- col
          else this_col <- altcol

          lines(xpos[indexes[[i]]], lod[indexes[[i]]],
                lwd=lwd, col=this_col, ...)
        }

        # add box just in case
        box()

      }

    # make the plot
    plot_scan2_internal(snp_lod = snp_lod, snpinfo = snpinfo, max_peak = max_peak,lod=lod,  haplo_peak_drop = haplo_peak_drop,
                        map=map, add=add, gap=gap, ...)
  }


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

# boundaries of chromosomes in plot_scan1
# first row: left edges
# second row: right edges
map_to_boundaries <-
  function(map, gap)
  {
    if(length(map)==1)
      return(cbind(range(map[[1]], na.rm=TRUE)))

    # range of each chromosome
    chr_range <- lapply(map, range, na.rm=TRUE)

    # corresponding xpos, as matrix with two rows
    startend <- matrix(map_to_xpos(chr_range, gap), nrow=2)

    startend[1,] <- startend[1,] - gap/2
    startend[2,] <- startend[2,] + gap/2

    startend
  }

# convert map to list of indexes to LOD vector
map_to_index <-
  function(map)
  {
    if(length(map)==1) {
      map[[1]] <- seq(along=map[[1]])
      return(map)
    }

    lengths <- vapply(map, length, 0)
    split(seq_len(sum(lengths)), rep(seq(along=map), lengths))
  }
