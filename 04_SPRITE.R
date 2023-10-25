##### Convert SPRITE clusters to the sparse matrix format of FitHiC2 #####

# Full matrix to sparse matrix; modified from HiCcompare
full2sparse <- function(mat) {
  require(data.table)

  if (is.null(colnames(mat))) {
    stop("Please set the column names of the matrix to the start location for
         each bin. See ?full2sparse")
  }
  regions <- colnames(mat)
  sparse <- as.data.table(expand.grid(regions, regions))
  up.triangle <- c(upper.tri(mat, diag = TRUE))
  sparse <- cbind(sparse, c(mat), up.triangle)
  colnames(sparse) <- c("region1", "region2", "IF", "up.triangle")
  sparse <- subset(sparse, up.triangle == TRUE & IF != 0)
  sparse <- sparse[, 1:3, with = FALSE]
  return(sparse)
}


sparse2fithic <- function(clusterFile, binFile) {
  clusters <- readr::read_tsv(clusterFile, col_names = FALSE)
  bins <- readr::read_tsv(binFile, col_names = FALSE)
  colnames(clusters) <- as.list(bins)$X1

  fithic <- full2sparse(as.matrix(clusters))

  return(fithic)
}


# Human
m1 <- sparse2fithic(
  "data/interchrom_contacts/SPRITE_1Mb/human_inter_1Mb_none_1000_raw.txt",
  "data/interchrom_contacts/SPRITE_1Mb/hg19_1Mb_bins.txt"
)

readr::write_tsv(m1, "GM12878.1Mb.Rawobserved.SPRITE.FitHiC2", col_names = FALSE)


# Mouse
m2 <- sparse2fithic(
  "data/interchrom_contacts/SPRITE_1Mb/mouse_inter_1Mb_none_1000_raw.txt",
  "data/interchrom_contacts/SPRITE_1Mb/mm9_1Mb_bins.txt"
)

readr::write_tsv(m2, "mESC.1Mb.Rawobserved.SPRITE.FitHiC2", col_names = FALSE)



##### Get colocalization status using SPRITE data #####
library(GenomicInteractions)
load("Rdata/3-simulations.Rdata")

# Human SPRITE
human_SPRITE_1Mb <- makeGenomicInteractionsFromFile(
  "data/interchrom_contacts/SPRITE_1Mb/GM12878_SPRITE_1Mb_q0.05.bedpe",
  type = "bedpe",
  experiment_name = "hs sprite",
  description = "hs sprite"
)

x <- colocalize(
  "data/interchrom_contacts/SPRITE_1Mb/GM12878_SPRITE_1Mb_q0.05.bedpe",
  human_retro
)

human_retro[, "GM12878_SPRITE_1Mb"] <- "Noncolocalized"
human_retro[x, ][, "GM12878_SPRITE_1Mb"] <- "Colocalized"


# Human HiC
human_HiC_1Mb <- makeGenomicInteractionsFromFile(
  "data/interchrom_contacts/SPRITE_1Mb/GM12878_HiC_1Mb_q0.05.bedpe",
  type = "bedpe",
  experiment_name = "hs hic",
  description = "hs hic"
)

x <- colocalize(
  "data/interchrom_contacts/SPRITE_1Mb/GM12878_HiC_1Mb_q0.05.bedpe",
  human_retro
)

human_retro[, "GM12878_HiC_1Mb"] <- "Noncolocalized"
human_retro[x, ][, "GM12878_HiC_1Mb"] <- "Colocalized"


# Mouse SPRITE
mouse_SPRITE_1Mb <- makeGenomicInteractionsFromFile(
  "data/interchrom_contacts/SPRITE_1Mb/mESC_SPRITE_1Mb_q0.05.bedpe",
  type = "bedpe",
  experiment_name = "mm sprite",
  description = "mm sprite"
)


x <- colocalize(
  "data/interchrom_contacts/SPRITE_1Mb/mESC_SPRITE_1Mb_q0.05.bedpe",
  mouse_retro
)

mouse_retro[, "mESC_SPRITE_1Mb"] <- "Noncolocalized"
mouse_retro[x, ][, "mESC_SPRITE_1Mb"] <- "Colocalized"


# Mouse HiC
mouse_HiC_1Mb <- makeGenomicInteractionsFromFile(
  "data/interchrom_contacts/SPRITE_1Mb/mESC_HiC_1Mb_q.05.bedpe",
  type = "bedpe",
  experiment_name = "mm hic",
  description = "mm hic"
)


x <- colocalize(
  "data/interchrom_contacts/SPRITE_1Mb/mESC_HiC_1Mb_q.05.bedpe",
  mouse_retro
)

mouse_retro[, "mESC_HiC_1Mb"] <- "Noncolocalized"
mouse_retro[x, ][, "mESC_HiC_1Mb"] <- "Colocalized"


save.image("Rdata/4-SPRITE.Rdata")


##### Venn plot #####

library(VennDiagram)
library(tidyverse)

venn_plot <- function(a1, a2, ca) {
  grid.newpage()

  draw.pairwise.venn(
    area1 = a1,
    area2 = a2,
    cross.area = ca,
    category = c("HiC", "SPRITE"),
    fill = c("#FF000080", "#0000FF80"),
    lty = "blank",
    cex = 2,
    cat.cex = 2,
    ext.pos = 180,
    cat.dist = c(0, 0),
  )
}


a1 <- nrow(human_retro %>% filter(GM12878_HiC_1Mb == "Colocalized"))
a2 <- nrow(human_retro %>% filter(GM12878_SPRITE_1Mb == "Colocalized"))
ca <- nrow(human_retro %>%
  filter(GM12878_HiC_1Mb == "Colocalized", GM12878_SPRITE_1Mb == "Colocalized"))

venn_plot(a1, a2, ca)

# Save plot as "Venn_plot_GM12878_1Mb.pdf", width = 9, height = 9


a1 <- nrow(mouse_retro %>% filter(mESC_HiC_1Mb == "Colocalized"))
a2 <- nrow(mouse_retro %>% filter(mESC_SPRITE_1Mb == "Colocalized"))
ca <- nrow(mouse_retro %>%
  filter(mESC_HiC_1Mb == "Colocalized", mESC_SPRITE_1Mb == "Colocalized"))

venn_plot(a1, a2, ca)

# Save plot as "Venn_plot_mESC_1Mb.pdf", width = 9, height = 9



##### Simulation #####

library(doParallel)
# Setup backend to use many processors
totalCores <- detectCores()

# Leave one core to avoid overload your computer
cluster <- makeCluster(totalCores[1] - 4)
registerDoParallel(cluster)

sim_hs <- run_sim(
  human_SPRITE_1Mb, human_retro, hg19_chr_sizes,
  human_genes, human_pcgenes, 1000
)

sim_mm <- run_sim(
  mouse_SPRITE_1Mb, mouse_retro, mm9_chr_sizes,
  mouse_genes, mouse_pcgenes, 1000
)


simTotal <- cbind(sim_hs, sim_mm)
names(simTotal) <- c("GM12878_SPRITE_1Mb", "mESC_SPRITE_1Mb")
simTotal[, "Simulation"] <- as.vector(sapply(1:5, function(i) rep(i, 1000)))

readr::write_tsv(simTotal, "results/simulation_SPRITE_1Mb.tsv")


##### Plot simulation results #####
source("simulation_code.R")
source("theme_publication.R")
library(tidyverse)
load("Rdata/4-SPRITE.Rdata")


sim_results <- read_tsv("results/simulation_SPRITE_1Mb.tsv")
sim_results$Simulation <- as.factor(sim_results$Simulation)


# GM12878
sim_plot(human_retro, "GM12878_SPRITE_1Mb", sim_results, 16, 1, 30, 1000)
ggsave("plots/figS4B.pdf", width = 8, height = 6)

# mESC
sim_plot(mouse_retro, "mESC_SPRITE_1Mb", sim_results, 16, 1, 30, 1000)
ggsave("plots/figS32B.pdf", width = 8, height = 6)
