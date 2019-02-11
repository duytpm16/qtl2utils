
plot_genes2 <-
  function(genes,  scan1output, snpinfo, snp_lod, max_peak, haplo_peak_drop = NULL,
           minrow=4, padding=0.2, text_cex = .2,
           colors=c("black", "red3", "green4", "blue3", "darkorange2"),
           scale_pos=1, start_field="start", stop_field="stop",
           strand_field="strand", name_field="Name", ...)
  {
    if(is.null(genes)) stop("genes is NULL")

    # make sure the columns are there
    fields <- c(start_field, stop_field, strand_field, name_field)
    fields_found <- fields %in% colnames(genes)
    if(!all(fields_found)) {
      stop("Columns not found: ", paste(fields[!fields_found], collapse=", "))
    }


    # grab just the start and stop
    start <- genes[,start_field]
    end <- genes[,stop_field]

    # drop genes with missing start or stop
    missing_pos <- is.na(start) | is.na(end)
    if(any(missing_pos)) {
      warning("Dropping ", sum(missing_pos), " rows with missing positions")
      genes <- genes[!missing_pos, , drop=FALSE]
    }

    # make sure genes are ordered by their start values
    if(any(diff(start) < 0))
      genes <- genes[order(start, stop),]

    # grab data
    start <- genes[,start_field]*scale_pos # convert to Mbp
    end <- genes[,stop_field]*scale_pos    # convert to Mbp
    strand <- as.character(genes[,strand_field])
    name <- as.character(genes[,name_field])

    plot_genes_internal <-
      function(snp_lod, haplo_peak_drop, max_peak, scan1output,
               snpinfo, xlab=NULL, xaxs="i",
               bgcolor="gray92", xat=NULL,
               mgp=c(1.6,0.2,0), xlim=NULL,
               vlines=NULL, vlines_col="white",
               vlines_lwd=1, vlines_lty=1,
               xaxt="s", ylab="", ...)
      {
        if(is.null(xlab)) {
          if(length(unique(genes$chr)) == 1)
            xlab <- paste("Chr", genes$chr[1], "position (Mbp)")
          else xlab <- "Position (Mbp)"
        }

        if(is.null(xlim)) {
          xlim <- range(c(start, end), na.rm=TRUE)
        }


        plot(0, 0, type="n",
             xlim=xlim,   xlab=xlab, xaxs=xaxs, xaxt="n",
             ylim=c(1,0), ylab=ylab, yaxs="i", yaxt="n", mgp=mgp, ...)


        max_snps <- find_max_snps(snp_pos = snpinfo$pos, snp_id = snpinfo$snp_id, snp_lod = snp_lod, max_peak = max_peak, haplo_peak_drop = NULL)
        snp_x <- c(max_snps$pos)


        # gray background
        u <- par("usr")
        rect(u[1], u[3], u[2], u[4], col=bgcolor, border="black")

        # axis
        if(is.null(xat)) xat <- pretty(xlim)
        if((length(xat) > 1 || !is.na(xat)) && xaxt != "n")
          axis(side=1, at=xat, mgp=mgp, tick=FALSE)

        # vertical lines
        if(is.null(vlines)) vlines <- xat
        if(length(vlines) > 1 || !is.na(vlines))
          abline(v=vlines, col=vlines_col, lwd=vlines_lwd, lty=vlines_lty)
        abline(v=max_peak)
        if(!is.null(haplo_peak_drop)){
          abline(v = haplo_peak_drop[1], lty = 3)
          abline(v = haplo_peak_drop[2], lty = 3)
        }
        abline(v=snp_x, col = "violetred")
        box()
      }

    plot_genes_internal(snp_lod = snp_lod, haplo_peak_drop = haplo_peak_drop,
                        max_peak = max_peak, scan1output = scan1output, snpinfo = snpinfo,...)

    # drop genes that are not in plotting region
    u <- par("usr")
    omit <- (end < u[1] | start > u[2])
    if(any(omit)) {
      keep <- !omit
      start <- start[keep]
      end <- end[keep]
      strand <- strand[keep]
      name <- name[keep]
    }
    if(length(start) == 0) # no genes to plot
      return(invisible(NULL))

    # missing names: use ?
    name[is.na(name)] <- "?"

    # arrow annotation re direction, to place after gene name
    dir_symbol <- rep(' ', length(name))
    right <- !is.na(strand) & strand == "+"
    if(any(right))
      dir_symbol[right] <- expression(phantom('') %->% phantom(''))
    left <- !is.na(strand) & strand == "-"
    if(any(left))
      dir_symbol[left] <- expression(phantom('') %<-% phantom(''))

    # initial determination of text size
    maxy <- minrow
    height <- 1/maxy
    text_cex <- text_cex

    # adjust text size and determine vertical location of genes
    for(it in 1:2) { # go through all of this twice

      while(max(abs(strheight(name, cex=text_cex))) > height*(1-padding)) {
        text_cex <- text_cex * 0.99
      }

      # horizontal padding
      space <- strwidth(' ', cex=text_cex)

      # figure out how to arrange genes vertically
      #     + number of rows of genes
      # (function defined in src/arrange_genes.cpp)
      y <- arrange_genes(start, end + space + strwidth(name, cex=text_cex) + strwidth(dir_symbol, cex=text_cex))

      maxy <- max(c(y, minrow))
      height <- 1/maxy
    }

    ypos <- seq(height/2, by=height, length=maxy)
    y <- ypos[y]
    rect_height <- height*(1-padding)
    rect_top <- y - rect_height/2
    rect_bottom <- y + rect_height/2

    colors <- rep(colors, length(y))

    for(i in seq(along=start)) {
      rect(start[i], rect_top[i],
           end[i],   rect_bottom[i],
           col=colors[i], border=colors[i],
           lend=1, ljoin=1)
      text(end[i] + space, y[i],
           name[i], adj=c(0, 0.5), col=colors[i],
           cex=text_cex)
      if(!is.na(strand[i]) && (strand[i] == "+" || strand[i] == '-'))
        text(end[i] + space + strwidth(name[i], cex=text_cex), y[i],
             dir_symbol[i], adj=c(0, 0.5), col=colors[i],
             cex=text_cex)
    }
  }


