##### mESC seqFISH+ data (200 kb) #####

library(GenomicInteractions)
library(tidyverse)

mouse_retro_mm10 <- read_tsv("data/interchrom_contacts/FISH/mouse_retro_inter_mm10.bedpe") %>%
  dplyr::rename(retroChr = `#retroChr`)


# Retro-parent pairs
mmRetroParent <- makeGenomicInteractionsFromFile(
  "data/interchrom_contacts/FISH/mouse_retro_inter_mm10.bedpe",
  type = "bedpe",
  experiment_name = "mm retro",
  description = "mm retro pairs"
)


# HiC interactions
mESC_HiC <- makeGenomicInteractionsFromFile(
  "data/interchrom_contacts/FISH/mESC_HiC_200kb_top5M.bedpe",
  type = "bedpe",
  experiment_name = "HiC 200kb",
  description = "mESC"
)

colocalPairs <- findOverlaps(mmRetroParent, mESC_HiC)
uniq <- unique(queryHits(colocalPairs))

mouse_retro_mm10$mESC_HiC <- "Noncolocalized"
mouse_retro_mm10[uniq, ]$mESC_HiC <- "Colocalized"

table(mouse_retro_mm10$mESC_HiC) # 1465


# seqFISH+ distances
seqFISH <- makeGenomicInteractionsFromFile(
  "data/interchrom_contacts/FISH/mESC_seqFISH+_200kb_interchr_distance_mat.bedpe",
  type = "bedpe",
  experiment_name = "seqfish 200kb",
  description = "mESC"
)

m <- mergeByOverlaps(mmRetroParent, seqFISH, ignore.strand = TRUE)
df <- as.data.frame(m) %>% select(mmRetroParent.NA., seqFISH.counts)
names(df) <- c("retroID", "distance")
df <- df %>%
  group_by(retroID) %>%
  summarise(dist_um = mean(distance))

mouse_retro_mm10 <- mouse_retro_mm10 %>% left_join(df, by = "retroID")

# Plot
source("theme_publication.R")
source("format_pvalues.R")
library(ggsignif)

p1 <- ggplot(mouse_retro_mm10, aes(x = mESC_HiC, y = dist_um)) +
  geom_violin(aes(color = mESC_HiC, fill = mESC_HiC), scale = 1.6, alpha = 1.0) +
  geom_boxplot(width = 0.2, outlier.shape = NA, size = 1.5) +
  geom_signif(
    comparisons = list(
      c("Colocalized", "Noncolocalized")
    ),
    map_signif_level = pform, size = 0.8, textsize = 8
  ) +
  labs(x = "", y = "Spatial distance (μm)", title = "seqFISH+") +
  theme_publication(baseSize = 24) +
  theme(
    legend.position = "none",
    legend.title = element_blank()
  )



########## Simulations for top 10% seqFISH+ in mESC ##########

load("3-simulations.Rdata")

library(tidyverse)
mouse_retro_mm10 <- read_tsv("data/interchrom_contacts/FISH/mouse_retro_inter_mm10.bedpe")

library(GenomicInteractions)
mESC_seqFISH <- makeGenomicInteractionsFromFile(
  "data/interchrom_contacts/FISH/mESC_seqFISH+_200kb_top0.1.bedpe",
  type = "bedpe",
  experiment_name = "seqfish 200kb",
  description = "mESC"
)

mm10_chr_sizes <- read_tsv("data/interchrom_contacts/FISH/mm10.chrom.sizes",
  col_names = F
)
colnames(mm10_chr_sizes) <- c("chrom", "size")


mouse_genes_mm10 <- get_gene_coords(
  "data/annotations/Mus_musculus.GRCm38.102.gtf.gz",
  mouse_chroms
)
mouse_pcgenes_mm10 <- get_pcgene_coords(
  "data/annotations/Mus_musculus.GRCm38.102.gtf.gz",
  mouse_chroms
)


# get co-localization status
x <- colocalize(
  "data/interchrom_contacts/FISH/mESC_seqFISH+_200kb_top0.1.bedpe",
  mouse_retro_mm10
)

mouse_retro_mm10[, "mESC_seqFISH"] <- "Noncolocalized"
mouse_retro_mm10[x, ][, "mESC_seqFISH"] <- "Colocalized"

save.image("Rdata/8-mESC_seqFISH.Rdata")


# run simulations

library(doParallel)
# Setup backend to use many processors
totalCores <- detectCores()

# Leave one core to avoid overload your computer
cluster <- makeCluster(totalCores[1] - 4)
registerDoParallel(cluster)


n <- 1000

sim_seqfish <- run_sim(
  mESC_seqFISH, mouse_retro_mm10, mm10_chr_sizes,
  mouse_genes, mouse_pcgenes, n
)

# Save simulation results
names(sim_seqfish) <- "mESC_seqFISH"

sim_seqfish[, "Simulation"] <- as.vector(sapply(1:5, function(i) rep(i, n)))
readr::write_tsv(sim_seqfish, "results/simulation_mESC_seqFISH.tsv")

### Plot simulation results
load("Rdata/8-mESC_seqFISH.Rdata")
sim_seqfish <- read_tsv("results/simulation_mESC_seqFISH.tsv")
sim_seqfish$Simulation <- as.factor(sim_seqfish$Simulation)

p2 <- sim_plot(mouse_retro_mm10, "mESC_seqFISH", sim_seqfish, 24, 2, 25, 1000) +
  theme(
    plot.title = element_blank(),
    legend.position = "top"
  )


# Save plot
library(ggpubr)
ggarrange(p1, p2,
  labels = "AUTO",
  font.label = list(size = 24, color = "black", face = "bold", family = "sans")
)

ggsave("plots/figS33.pdf", width = 20, height = 10, device = cairo_pdf)
