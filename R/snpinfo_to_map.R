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
