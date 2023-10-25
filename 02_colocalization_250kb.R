########## Colocalization status at 250 kb for all species ##########

load("Rdata/1-retrocopy_info.Rdata")


colocalize <- function(HiC, retroInfo) {
  require(GenomicInteractions)

  GI <- makeGenomicInteractionsFromFile(
    HiC,
    type = "bedpe",
    experiment_name = "HiC 250kb",
    description = HiC
  )

  R <- GRanges(
    seqnames = retroInfo$retroChr,
    ranges = IRanges(
      start = retroInfo$retroStart,
      end = retroInfo$retroEnd
    ),
    strand = retroInfo$retroStrand,
    names = retroInfo$retroID
  )

  P <- GRanges(
    seqnames = retroInfo$parentChr,
    ranges = IRanges(
      start = retroInfo$parentStart,
      end = retroInfo$parentEnd
    ),
    strand = retroInfo$parentStrand,
    names = retroInfo$parentID
  )

  RP <- GenomicInteractions(R, P)

  colocalPairs <- findOverlaps(RP, GI)
  uniq <- unique(queryHits(colocalPairs))

  return(uniq)
}



##### Human cell lines (hg19) #####
cellLines <- c("GM12878", "HMEC", "HUVEC", "IMR90", "K562", "KBM7", "NHEK")

for (cellLine in cellLines) {
  x <- colocalize(paste0(
    "data/interchrom_contacts/HiC_250kb/",
    cellLine, "_HiC_250kb_top5M.bedpe"
  ), human_retro)

  human_retro[, cellLine] <- "Noncolocalized"
  human_retro[x, ][, cellLine] <- "Colocalized"
}


##### Mouse CH12-LX cell line (mm9) #####
x <- colocalize(
  "data/interchrom_contacts/HiC_250kb/CH12-LX_HiC_250kb_top5M.bedpe",
  mouse_retro
)

mouse_retro[, "CH12.LX"] <- "Noncolocalized"
mouse_retro[x, ][, "CH12.LX"] <- "Colocalized"


##### Chimpanzee iPSC (panTro6) #####
x <- colocalize(
  "data/interchrom_contacts/HiC_250kb/Chimp_HiC_250kb_top5M.bedpe",
  chimp_retro
)

chimp_retro[, "Chimp"] <- "Noncolocalized"
chimp_retro[x, ][, "Chimp"] <- "Colocalized"


##### Marmoset iPSC (calJac3) #####

x <- colocalize(
  "data/interchrom_contacts/HiC_250kb/Marmoset_HiC_250kb_top5M.bedpe",
  marmoset_retro
)

marmoset_retro[, "Marmoset"] <- "Noncolocalized"
marmoset_retro[x, ][, "Marmoset"] <- "Colocalized"


##### Macaque CP and GZ (rheMac8) #####

for (s in c("Macaque_CP", "Macaque_GZ")) {
  x <- colocalize(
    paste0(
      "data/interchrom_contacts/HiC_250kb/Macaque_", s,
      "_HiC_250kb_top5M.bedpe"
    ),
    macaque_retro
  )

  macaque_retro[, s] <- "Noncolocalized"
  macaque_retro[x, ][, s] <- "Colocalized"
}


##### Cow ear skin (UMD3.1) #####

x <- colocalize(
  "data/interchrom_contacts/HiC_250kb/Cow_HiC_250kb_top5M.bedpe",
  cow_retro
)

cow_retro[, "Cow"] <- "Noncolocalized"
cow_retro[x, ][, "Cow"] <- "Colocalized"


##### Dog ear skin (CanFam3.1) #####

x <- colocalize(
  "data/interchrom_contacts/HiC_250kb/Dog_HiC_250kb_top5M.bedpe",
  dog_retro
)

dog_retro[, "Dog"] <- "Noncolocalized"
dog_retro[x, ][, "Dog"] <- "Colocalized"



# Save data
save.image("Rdata/2-colocalization_250kb.Rdata")
