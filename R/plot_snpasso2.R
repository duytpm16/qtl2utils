#' Plot SNP association#'
plot_snpasso2 <-
  function(scan1output, snpinfo, max_peak, haplo_peak_drop = NULL, haplo_scan = NULL, genes=NULL, lodcolumn=1, show_all_snps=TRUE, chr=NULL,
           add=FALSE, drop_hilit=NA, col_hilit="violetred", col="darkslateblue",
           gap=25, minlod=0, text_cex = .2, ...)
  {
    if(is.null(scan1output)) stop("scan1output is NULL")
    if(is.null(snpinfo)) stop("snpinfo is NULL")


    # pull out lod scores
    if(length(lodcolumn)==0) stop("lodcolumn has length 0")
    if(length(lodcolumn) > 1) { # If length > 1, take first value
      warning("lodcolumn should have length 1; only first element used.")
      lodcolumn <- lodcolumn[1]
    }
    if(is.character(lodcolumn)) { # turn column name into integer
      tmp <- match(lodcolumn, colnames(scan1output))
      if(is.na(tmp)) stop('lodcolumn "', lodcolumn, '" not found')
      lodcolumn <- tmp
    }
    if(lodcolumn < 1 || lodcolumn > ncol(scan1output))
      stop("lodcolumn [", lodcolumn, "] out of range (should be in 1, ..., ", ncol(scan1output), ")")
    scan1output <- scan1output[,lodcolumn,drop=FALSE]

    if(!is.null(chr)) { # subset by chr
      if(!any(chr %in% snpinfo$chr)) stop("None of the chr found in snpinfo")
      if(!all(chr %in% snpinfo$chr)) {
        warning("Some chr not found: ", paste(chr[!(chr %in% snpinfo$chr)], collapse=", "))
      }
      snpinfo <- snpinfo[snpinfo$chr %in% chr,,drop=FALSE]
      scan1output <- scan1output[rownames(scan1output) %in% snpinfo$snp,,drop=FALSE]
    }

    if(!is.null(genes)) {
      if(length(unique(snpinfo$chr)) > 1) {
        warning("genes ignored if there are multiple chromosomes")
      } else {
        return( plot_snpasso_and_genes2(scan1output, snpinfo, max_peak = max_peak, haplo_peak_drop = haplo_peak_drop, haplo_scan = haplo_scan,
                                        show_all_snps=show_all_snps,
                                        drop_hilit=drop_hilit, col_hilit=col_hilit,
                                        col=col, gap=gap, minlod=minlod, genes=genes, text_cex = text_cex, ...) )
      }
    }

    if(nrow(scan1output) == nrow(snpinfo) && all(rownames(scan1output) == snpinfo$snp)) {
      show_all_snps <- FALSE

      snpinfo <- snpinfo[scan1output[,1]>=minlod, , drop=FALSE]
      scan1output <- scan1output[scan1output[,1]>=minlod, 1, drop=FALSE]

      # skip the index business
      # snpinfo -> map
      map <- tapply(seq_len(nrow(snpinfo)), factor(snpinfo$chr, unique(snpinfo$chr)),
                    function(i) setNames(snpinfo$pos[i], snpinfo$snp[i]))
    }
    else {
      snpinfo_spl <- split(snpinfo, factor(snpinfo$chr, unique(snpinfo$chr)))

      uindex <- unlist(lapply(snpinfo_spl, function(a) unique(a$index)))
      if(length(uindex) != nrow(scan1output)) {
        stop("Something is wrong with snpinfo$index.\n",
             "      no. unique indexes [",
             length(uindex), "] != nrow(scan1output) [",
             nrow(scan1output), "].")
      }

      for(i in seq_along(snpinfo_spl)) {
        uindex <- unique(snpinfo_spl[[i]]$index)
        if(any(snpinfo_spl[[i]]$index[uindex] != uindex)) {
          stop("Something is wrong with snpinfo$index.\n",
               "      on each chr, index[u] should == u for each index u")
        }
      }

      map <- snpinfo_to_map(snpinfo)

      if(show_all_snps) {
        tmp <- expand_snp_results(scan1output, map, snpinfo)
        scan1output <- tmp$lod
        map <- tmp$map
      }
    }

    # maximum LOD
    maxlod <- max(unclass(scan1output)[,1], na.rm=TRUE)

    # set values < minlod to NA so they're not plotted
    scan1output[scan1output < minlod] <- NA

    # internal function to give defaults to hidden graphics parameters
    plot_snpasso_internal2 <-
      function(max_peak, haplo_peak_drop = NULL, haplo_scan = NULL, pch=16, cex=0.5, ylim=NULL, bgcolor="gray90",
               altbgcolor="gray85",
               drop_hilit=NA, col_hilit="violetred",
               drop.hilit=NULL, col.hilit=NULL,text_cex = text_cex, ...)
      {
        if(!is.null(drop.hilit)) {
          warning("drop.hilit is deprecated; use drop_hilit")
          drop_hilit <- drop.hilit
        }

        if(!is.null(col.hilit)) {
          warning("col.hilit is deprecated; use col_hilit")
          col_hilit <- col.hilit
        }

        if(is.null(ylim))
          ylim <- c(0, maxlod*1.02)

        if(!is.na(drop_hilit) && !is.null(drop_hilit)) {
          col <- c(col, col_hilit)[(scan1output >= maxlod-drop_hilit)+1]
        }

        plot_scan2(x = scan1output, snp_lod = scan1output, snpinfo = snpinfo, max_peak = max_peak, haplo_peak_drop = haplo_peak_drop,
                   map = map, lodcolumn=1, bgcolor=bgcolor, altbgcolor=altbgcolor, ylim=ylim,
                   gap=gap, add=add, col = col, type="p", cex=cex, pch=pch,...)
      }

    plot_snpasso_internal2( max_peak = max_peak, haplo_peak_drop = haplo_peak_drop, drop_hilit=drop_hilit, col_hilit=col_hilit,text_cex = text_cex, ...)
  }

# expand snp association results according to snpinfo
expand_snp_results <-
  function(snp_results, map, snpinfo)
  {
    snpinfo <- split(snpinfo, factor(snpinfo$chr, unique(snpinfo$chr)))

    if(length(map) != length(snpinfo))
      stop("length(map) [", length(map), "] != length(snpinfo) [",
           length(snpinfo), "]")

    if(nrow(snp_results) != length(unlist(map)))
      stop("nrow(snp_results) [", nrow(snp_results), "] != length(unlist(map)) [",
           length(unlist(map)), "]")

    cnames <- rep(names(map), vapply(map, length, 0))
    lodindex <- split(seq_len(nrow(snp_results)), factor(cnames, unique(cnames)))

    result <- NULL
    for(i in seq(along=map)) {
      revindex <- rev_snp_index(snpinfo[[i]])

      map[[i]] <- snpinfo[[i]]$pos
      names(map[[i]]) <- snpinfo[[i]]$snp
      this_result <- unclass(snp_results)[lodindex[[i]],,drop=FALSE][revindex,,drop=FALSE]
      rownames(this_result) <- snpinfo[[i]]$snp

      result <- rbind(result, this_result)
    }

    list(lod=result,
         map=map)
  }


# snpinfo to map
snpinfo_to_map <-
  function(snpinfo)
  {
    uindex <- sort(unique(snpinfo$index))
    if(any(snpinfo$index < 1 | snpinfo$index > nrow(snpinfo)))
      stop("snpinfo$index values outside of range [1, ",
           nrow(snpinfo), "]")

    uchr <- unique(snpinfo$chr)
    chr <- factor(snpinfo$chr, levels=uchr)

    map <- split(snpinfo$pos, chr)
    snp <- split(snpinfo$snp, chr)
    index <- split(snpinfo$index, chr)
    for(i in seq_along(map)) {
      u <- unique(index[[i]])
      map[[i]] <- map[[i]][u]
      names(map[[i]]) <- snp[[i]][u]
    }

    names(map) <- uchr

    map
  }

# reverse index
rev_snp_index <-
  function(snpinfo)
  {
    index_spl <- split(seq_len(nrow(snpinfo)), snpinfo$index)
    revindex <- rep(seq_along(index_spl), vapply(index_spl, length, 1))
    revindex[unlist(index_spl)] <- revindex

    revindex
  }

