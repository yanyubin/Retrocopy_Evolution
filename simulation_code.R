########## Simulation code ##########

# shuffleBed
shuffleBed <- function(gr) {
  require(rtracklayer)
  require(GenomicInteractions)

  stopifnot(class(gr) == "GRanges")
  if (nchar(Sys.which("bedtools")) == 0) stop("shuffleBed must be in the PATH")
  if (any(is.na(seqlengths(gr)))) stop("seqlengths of grange must be defined")

  genome_file <- tempfile()
  writeLines(paste(names(seqlengths(gr)), seqlengths(gr), sep = "\t"), genome_file)

  bed_file <- tempfile()
  writeLines(paste(as.character(seqnames(gr)), start(gr) - 1, end(gr),
    sep = "\t",
    strand(gr)
  ), bed_file)

  shuffle_file <- tempfile()
  system(paste(
    "bedtools shuffle -i", bed_file, "-g", genome_file,
    "-chrom", ">", shuffle_file
  ))

  shuffled <- import.bed(shuffle_file)
  seqlevels(shuffled) <- seqlevels(gr)
  seqlengths(shuffled) <- seqlengths(gr)

  unlink(c(genome_file, bed_file, shuffle_file))

  return(shuffled)
}


# Simulation 1: random pairs with matched chromosome and length
simulation1 <- function(SI, retroInfo, chromSizes, n) {
  require(GenomicInteractions)
  shuffleBed <- shuffleBed

  R <- df2grange(retroInfo, "retrocopy")
  seqlengths(R) <- chromSizes$size[match(
    names(seqlengths(R)),
    chromSizes$chrom
  )]

  P <- df2grange(retroInfo, "parent")
  seqlengths(P) <- chromSizes$size[match(
    names(seqlengths(P)),
    chromSizes$chrom
  )]

  foreach(i = 1:n, .combine = rbind) %dopar% {
    shuffleR <- shuffleBed(R)
    shuffleP <- shuffleBed(P)
    shuffleRP <- GenomicInteractions(shuffleR, shuffleP, counts = 1)
    paired_hits <- findOverlaps(shuffleRP, SI)
    r <- length(unique(queryHits(paired_hits)))

    return(r)
  }
}

# simulation1(GM12878_HiC, human_retro, hg19_chr_sizes, 1000)


# Simulation 2: retrocopy-random gene pairs, including non-coding genes
simulation2 <- function(SI, retroInfo, genes, n) {
  require(GenomicInteractions)
  df2grange <- df2grange

  R <- df2grange(retroInfo, "retrocopy")

  dtR <- data.table::data.table(retroInfo)
  dtG <- data.table::data.table(genes)

  foreach(i = 1:n, .combine = rbind) %dopar% {
    # randomly sample one gene with matching chromosome as the parental gene for each retro
    dtP <- dtG[dtR,
      on = .(parentChr),
      {
        ri <- sample(.N, 1L)
        .(
          parentStart = parentStart[ri], parentEnd = parentEnd[ri],
          parentStrand = parentStrand[ri], parentID = parentID[ri]
        )
      },
      by = .EACHI
    ]

    parentTmp <- df2grange(dtP, "parent")

    shuffleRP <- GenomicInteractions(R, parentTmp, counts = 1)
    paired_hits <- findOverlaps(shuffleRP, SI)
    r <- length(unique(queryHits(paired_hits)))

    return(r)
  }
}

# simulation2(GM12878_HiC, human_retro, human_genes, 1000)


# Simulation 3: parental gene-random "retrocopy" pairs
simulation3 <- function(SI, retroInfo, chromSizes, n) {
  require(GenomicInteractions)
  shuffleBed <- shuffleBed
  df2grange <- df2grange

  R <- df2grange(retroInfo, "retrocopy")
  seqlengths(R) <- chromSizes$size[match(
    names(seqlengths(R)),
    chromSizes$chrom
  )]

  P <- df2grange(retroInfo, "parent")

  foreach(i = 1:n, .combine = rbind) %dopar% {
    shuffleR <- shuffleBed(R)
    shuffleRP <- GenomicInteractions(shuffleR, P, counts = 1)
    paired_hits <- findOverlaps(shuffleRP, SI)
    r <- length(unique(queryHits(paired_hits)))

    return(r)
  }
}

# simulation3(GM12878_HiC, human_retro, hg19_chr_sizes, 1000)


# Simulation 4: retrocopy-random gene pairs, protein-coding gene only
simulation4 <- function(SI, retroInfo, pcgenes, n) {
  require(GenomicInteractions)
  df2grange <- df2grange

  R <- df2grange(retroInfo, "retrocopy")

  dtR <- data.table::data.table(retroInfo)
  dtG <- data.table::data.table(pcgenes)

  foreach(i = 1:n, .combine = rbind) %dopar% {
    # randomly sample one gene with matching chromosome as the parental gene for each retro
    dtP <- dtG[dtR,
      on = .(parentChr),
      {
        ri <- sample(.N, 1L)
        .(
          parentStart = parentStart[ri], parentEnd = parentEnd[ri],
          parentStrand = parentStrand[ri], parentID = parentID[ri]
        )
      },
      by = .EACHI
    ]

    parentTmp <- df2grange(dtP, "parent")

    shuffleRP <- GenomicInteractions(R, parentTmp, counts = 1)
    paired_hits <- findOverlaps(shuffleRP, SI)
    r <- length(unique(queryHits(paired_hits)))

    return(r)
  }
}

# simulation4(GM12878_HiC, human_retro, human_pcgenes, 1000)


# Simulation 5: shuffle retrocopy-parental gene pairs
simulation5 <- function(SI, retroInfo, n) {
  require(GenomicInteractions)
  df2grange <- df2grange

  R <- df2grange(retroInfo, "retrocopy")

  foreach(i = 1:n, .combine = rbind) %dopar% {
    Temp <- retroInfo[sample(1:nrow(retroInfo)), ]

    parentTemp <- df2grange(Temp, "parent")

    shuffleRP <- GenomicInteractions(R, parentTemp, counts = 1)
    paired_hits <- findOverlaps(shuffleRP, SI)
    r <- length(unique(queryHits(paired_hits)))

    return(r)
  }
}

# simulation5(GM12878_HiC, human_retro, 1000)


# a function to plot simulation results
sim_plot <- function(retroInfo, cellLine, sim, baseSize, arrowSize, arrowLen, nSim) {
  require(tidyverse)

  colNum <- retroInfo %>%
    filter(eval(parse(text = cellLine)) == "Colocalized") %>%
    nrow()
  pval <- sim %>%
    filter(Simulation == "5", eval(parse(text = cellLine)) >= colNum) %>%
    nrow()
  pval <- if (pval == 0) {
    bquote(italic("p") * " <" ~ .(1 / nSim))
  } else {
    bquote(italic("p") * " =" ~ .(pval / nSim))
  }
  n_retro <- nrow(retroInfo)

  ggplot(data = sim, aes(x = eval(parse(text = cellLine)) / n_retro, fill = Simulation)) +
    geom_density(alpha = .5) +
    geom_segment(
      aes(
        x = colNum / n_retro, xend = colNum / n_retro,
        y = arrowLen, yend = 0, color = "red"
      ),
      size = arrowSize, show.legend = FALSE,
      arrow = arrow(length = unit(0.2, "cm"))
    ) +
    labs(
      title = cellLine, subtitle = pval,
      x = "Colocalization frequency", y = "Density",
      fill = "Simulation"
    ) +
    theme_publication(baseSize = baseSize) +
    theme(legend.position = "bottom")
}
