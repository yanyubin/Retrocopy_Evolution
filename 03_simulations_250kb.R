########## Simulations at 250 kb ##########

##### Generate dataframes with protein coding genes #####

# A function to get protein-coding gene coordinates from gtf file
get_pcgene_coords <- function(gtf, chroms) {
  require(ensembldb)
  require(dplyr)

  db <- ensDbFromGtf(gtf = gtf)
  edb <- EnsDb(db)
  df <- as.data.frame(genes(edb, filter = GeneBiotypeFilter("protein_coding"))) %>%
    dplyr::select(seqnames:end, strand:gene_id) %>%
    dplyr::mutate(seqnames = paste0("chr", seqnames, sep = "")) %>%
    dplyr::filter(seqnames %in% chroms)

  colnames(df) <- c(
    "parentChr", "parentStart", "parentEnd",
    "parentStrand", "parentID"
  )

  return(df)
}


load("Rdata/2-colocalization_250kb.Rdata")

# Get protein-coding gene coordinates of 7 species
human_pcgenes <- get_pcgene_coords(
  "data/annotations/Homo_sapiens.GRCh37.65.gtf.gz",
  human_chroms
)

mouse_pcgenes <- get_pcgene_coords(
  "data/annotations/Mus_musculus.NCBIM37.65.gtf.gz",
  mouse_chroms
)

macaque_genes <- get_gene_coords(
  "data/annotations/Macaca_mulatta.Mmul_8.0.1.97.gtf.gz",
  macaque_chroms
)

macaque_pcgenes <- get_pcgene_coords(
  "data/annotations/Macaca_mulatta.Mmul_8.0.1.97.gtf.gz",
  macaque_chroms
)

marmoset_pcgenes <- get_pcgene_coords(
  "data/annotations/Callithrix_jacchus.C_jacchus3.2.1.80.gtf.gz",
  marmoset_chroms
)

cow_pcgenes <- get_pcgene_coords(
  "data/annotations/Bos_taurus.UMD3.1.94.gtf.gz",
  cow_chroms
)

dog_pcgenes <- get_pcgene_coords(
  "data/annotations/Canis_lupus_familiaris.CanFam3.1.104.gtf.gz",
  dog_chroms
)

# Read manually curated pc genes of chimpanzee

chimp_genes <- read.table("data/annotations/panTro6.ncbiRefSeq.genes.tsv")

names(chimp_genes) <- c(
  "parentChr", "parentStart", "parentEnd",
  "parentStrand", "parentID"
)
library(tidyverse)
chimp_genes <- chimp_genes %>%
  group_by(parentID) %>%
  summarise(
    parentChr = min(parentChr), parentStart = min(parentStart),
    parentEnd = max(parentEnd), parentStrand = min(parentStrand)
  )



##### Read chromosome sizes as dataframe #####
spp <- c("hg19", "mm9", "panTro6", "rheMac8", "calJac3", "Cow.UMD3.1", "CanFam3.1")

for (s in spp) {
  assign(
    paste0(s, "_chr_sizes"),
    read.table(paste0("data/refseqs/chrom_sizes/", s, ".chrom.sizes"),
      col.names = c("chrom", "size")
    )
  )
}



##### Read HiC interaction data #####
hic_data <- c(
  "GM12878", "HMEC", "HUVEC", "IMR90", "K562", "KBM7", "NHEK",
  "CH12.LX", "Chimp", "Macaque_CP", "Marmoset", "Cow", "Dog"
)

for (s in hic_data) {
  assign(
    paste0(s, "_HiC"),
    makeGenomicInteractionsFromFile(
      paste0(
        "data/interchrom_contacts/HiC_250kb/",
        s, "_HiC_250kb_top5M.bedpe"
      ),
      type = "bedpe",
      experiment_name = "HiC 250kb",
      description = s
    )
  )
}



##### Generate Grange objects from retrocopy dataframes #####

df2grange <- function(retroInfo, type = "retrocopy") {
  require(GenomicRanges)

  if (type == "retrocopy") {
    g <- GRanges(
      seqnames = retroInfo$retroChr,
      ranges = IRanges(
        start = retroInfo$retroStart,
        end = retroInfo$retroEnd
      ),
      strand = retroInfo$retroStrand,
      names = retroInfo$retroID
    )
  } else if (type == "parent") {
    g <- GRanges(
      seqnames = retroInfo$parentChr,
      ranges = IRanges(
        start = retroInfo$parentStart,
        end = retroInfo$parentEnd
      ),
      strand = retroInfo$parentStrand,
      names = retroInfo$parentID
    )
  }

  return(g)
}


# Save data
save.image("Rdata/3-simulations.Rdata")



# A function to run simulations
run_sim <- function(SI, retroInfo, chromSizes, genes, pcgenes, n) {
  print("Running simulation 1")
  sim1 <- as.data.frame(simulation1(SI, retroInfo, chromSizes, n))

  print("Running simulation 2")
  sim2 <- as.data.frame(simulation2(SI, retroInfo, genes, n))

  print("Running simulation 3")
  sim3 <- as.data.frame(simulation3(SI, retroInfo, chromSizes, n))

  print("Running simulation 4")
  sim4 <- as.data.frame(simulation4(SI, retroInfo, pcgenes, n))

  print("Running simulation 5")
  sim5 <- as.data.frame(simulation5(SI, retroInfo, n))

  sim <- rbind(sim1, sim2, sim3, sim4, sim5)

  return(sim)
}


##### Run simulations #####
load("Rdata/3-simulations.Rdata")

library(doParallel)
# Setup backend to use many processors
totalCores <- detectCores()

# Leave one core to avoid overload your computer
cluster <- makeCluster(totalCores[1] - 1)
registerDoParallel(cluster)

n <- 1000

for (s in hic_data[1:7]) {
  assign(paste0("sim_", s), run_sim(
    get(paste0(s, "_HiC")), human_retro, hg19_chr_sizes,
    human_genes, human_pcgenes, n
  ))
}

sim_CH12.LX <- run_sim(
  `CH12-LX_HiC`, mouse_retro, mm9_chr_sizes,
  mouse_genes, mouse_pcgenes, n
)

sim_Chimp <- run_sim(
  Chimp_HiC, chimp_retro, panTro6_chr_sizes,
  chimp_genes, chimp_genes, n
)

sim_Macaque_CP <- run_sim(
  Macaque_CP_HiC, macaque_retro, rheMac8_chr_sizes,
  macaque_genes, macaque_pcgenes, n
)

sim_Marmoset <- run_sim(
  Marmoset_HiC, marmoset_retro, calJac3_chr_sizes,
  marmoset_genes, marmoset_pcgenes, n
)

sim_Cow <- run_sim(
  Cow_HiC, cow_retro, Cow.UMD3.1_chr_sizes,
  cow_genes, cow_pcgenes, n
)

sim_Dog <- run_sim(
  Dog_HiC, dog_retro, CanFam3.1_chr_sizes,
  dog_genes, dog_pcgenes, n
)


# Save simulation results
sim <- cbind(
  sim_GM12878, sim_HMEC, sim_HUVEC, sim_IMR90, sim_K562, sim_KBM7, sim_NHEK,
  sim_CH12.LX, sim_Chimp, sim_Macaque_CP, sim_Marmoset, sim_Cow, sim_Dog
)
names(sim) <- hic_data


sim[, "Simulation"] <- as.vector(sapply(1:5, function(i) rep(i, n)))
readr::write_tsv(sim, "results/simulation_250kb.tsv")


# Stop cluster
stopCluster(cluster)



##### Plot simulation results #####
source("simulation_code.R")
source("theme_publication.R")
library(tidyverse)
load("Rdata/3-simulations.Rdata")

sim_results <- read_tsv("results/simulation_250kb_results.tsv")
sim_results$Simulation <- as.factor(sim_results$Simulation)


library(ggpubr)

sim_plot(human_retro, "GM12878", sim_results, 28, 2, 25, 1000)
ggsave("./plots/fig1C.pdf", width = 10, height = 8)


p1 <- sim_plot(human_retro, "HMEC", sim_results, 16, 1, 25, 1000)
p2 <- sim_plot(human_retro, "HUVEC", sim_results, 16, 1, 25, 1000)
p3 <- sim_plot(human_retro, "IMR90", sim_results, 16, 1, 25, 1000)
p4 <- sim_plot(human_retro, "K562", sim_results, 16, 1, 25, 1000)
p5 <- sim_plot(human_retro, "KBM7", sim_results, 16, 1, 25, 1000)
p6 <- sim_plot(human_retro, "NHEK", sim_results, 16, 1, 25, 1000)

ggarrange(p1, p2, p3, p4, p5, p6,
  common.legend = TRUE, labels = "AUTO",
  font.label = list(color = "black", size = 18)
)
ggsave("./plots/figS3.pdf", width = 12, height = 8)


# Other mammals
p1 <- sim_plot(mouse_retro, "CH12.LX", sim_results, 16, 1, 25, 1000)
p2 <- sim_plot(chimp_retro, "Chimp", sim_results, 16, 1, 20, 1000)
p3 <- sim_plot(macaque_retro, "Macaque_CP", sim_results, 16, 1, 25, 1000)
p4 <- sim_plot(marmoset_retro, "Marmoset", sim_results, 16, 1, 20, 1000)
p5 <- sim_plot(cow_retro, "Cow", sim_results, 16, 1, 15, 1000)
p6 <- sim_plot(dog_retro, "Dog", sim_results, 16, 1, 15, 1000)

ggarrange(p1, p2, p3, p4, p5, p6,
          common.legend = TRUE, labels = "AUTO",
          font.label = list(color = "black", size = 18)
)
ggsave("./plots/figS31.pdf", width = 12, height = 8)
