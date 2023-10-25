##### Compare GM12878 100kb, 250kb, and 500kb #####
load("Rdata/1-retrocopy_info.Rdata")

library(GenomicInteractions)

# retro-parent pairs
hsRetroParent <- makeGenomicInteractionsFromFile(
  "data/interchrom_contacts/GM12878/human_retro.bedpe",
  type = "bedpe",
  experiment_name = "hs retro",
  description = "hs retro pairs"
)


# 100 kb
HiC100 <- makeGenomicInteractionsFromFile(
  "data/interchrom_contacts/GM12878/GM12878_HiC_100kb_q0.05.bedpe",
  type = "bedpe",
  experiment_name = "HiC 100kb",
  description = "GM12878 100kb"
)

colocalPairs <- findOverlaps(hsRetroParent, HiC100)
uniq <- unique(queryHits(colocalPairs))
human_retro[, "GM12878_100kb"] <- "Noncolocalized"
human_retro[uniq, ][, "GM12878_100kb"] <- "Colocalized"


# 250 kb
HiC250 <- makeGenomicInteractionsFromFile(
  "data/interchrom_contacts/GM12878/GM12878_HiC_250kb_q0.05.bedpe",
  type = "bedpe",
  experiment_name = "HiC 250kb",
  description = "GM12878 250kb"
)

colocalPairs <- findOverlaps(hsRetroParent, HiC250)
uniq <- unique(queryHits(colocalPairs))
human_retro[, "GM12878_250kb"] <- "Noncolocalized"
human_retro[uniq, ][, "GM12878_250kb"] <- "Colocalized"


# 500 kb
HiC500 <- makeGenomicInteractionsFromFile(
  "data/interchrom_contacts/GM12878/GM12878_HiC_500kb_q0.05.bedpe",
  type = "bedpe",
  experiment_name = "HiC 500kb",
  description = "GM12878 500kb"
)

colocalPairs <- findOverlaps(hsRetroParent, HiC500)
uniq <- unique(queryHits(colocalPairs))
human_retro[, "GM12878_500kb"] <- "Noncolocalized"
human_retro[uniq, ][, "GM12878_500kb"] <- "Colocalized"

human_retro_GM12878 <- human_retro

rm(colocalPairs, uniq)
save.image("Rdata/10-GM12878_colocalization.Rdata")


# Venn plot
load("Rdata/10-GM12878_colocalization.Rdata")

table(human_retro$GM12878_100kb)
table(human_retro$GM12878_250kb)
table(human_retro$GM12878_500kb)

library(VennDiagram)
library(tidyverse)

n12 <- human_retro %>%
  filter(GM12878_100kb == "Colocalized", GM12878_250kb == "Colocalized") %>%
  nrow()

n23 <- human_retro %>%
  filter(GM12878_250kb == "Colocalized", GM12878_500kb == "Colocalized") %>%
  nrow()

n13 <- human_retro %>%
  filter(GM12878_100kb == "Colocalized", GM12878_500kb == "Colocalized") %>%
  nrow()

n123 <- human_retro %>%
  filter(GM12878_100kb == "Colocalized", GM12878_250kb == "Colocalized", GM12878_500kb == "Colocalized") %>%
  nrow()


grid.newpage()

venn.plot <- draw.triple.venn(
  area1 = 244,
  area2 = 1039,
  area3 = 1495,
  n12 = n12,
  n23 = n23,
  n13 = n13,
  n123 = n123,
  category = c("100 kb", "250 kb", "500 kb"),
  fill = c("#FF000080", "#0000FF80", "greenyellow"),
  lty = "blank",
  cex = 2,
  cat.cex = 2,
)

# save image as "figS44.pdf", width = 7, height = 7


##### Simulations at 100 kb and 500 kb in GM12878 ######

load("3-simulations.Rdata")

# Read HiC interaction data
library(GenomicInteractions)

res <- c("100kb", "500kb")
for (r in res) {
  assign(
    paste0("GM12878_HiC_", r),
    makeGenomicInteractionsFromFile(
      paste0(
        "data/interchrom_contacts/GM12878/GM12878_HiC_",
        r, "_q0.05.bedpe"
      ),
      type = "bedpe",
      experiment_name = "HiC",
      description = r
    )
  )
}


for (r in res) {
  x <- colocalize(paste0(
    "data/interchrom_contacts/GM12878/GM12878_HiC_",
    r, "_q0.05.bedpe"
  ), human_retro)

  cellLine <- paste0("GM12878_", r)

  human_retro[, cellLine] <- "Noncolocalized"
  human_retro[x, ][, cellLine] <- "Colocalized"
}


# run simulations

library(doParallel)
# Setup backend to use many processors
totalCores <- detectCores()

# Leave one core to avoid overload your computer
cluster <- makeCluster(totalCores[1] - 4)
registerDoParallel(cluster)


n <- 1000

for (r in res) {
  assign(paste0("sim_GM12878_", r), run_sim(
    get(paste0("GM12878_HiC_", r)), human_retro,
    hg19_chr_sizes, human_genes, human_pcgenes, n
  ))
}

# Save simulation results
sim <- cbind(sim_GM12878_100kb, sim_GM12878_500kb)
names(sim) <- c("GM12878_100kb", "GM12878_500kb")

sim[, "Simulation"] <- as.vector(sapply(1:5, function(i) rep(i, n)))
readr::write_tsv(sim, "results/simulation_GM12878_100kb_500kb.tsv")


# Plot simulation results
load("Rdata/3-simulations.Rdata")

sim <- read_tsv("results/simulation_GM12878_100kb_500kb.tsv")
sim$Simulation <- as.factor(sim$Simulation)

p1 <- sim_plot(human_retro_GM12878, "GM12878_100kb", sim, 24, 2, 50, 1000)
p2 <- sim_plot(human_retro_GM12878, "GM12878_500kb", sim, 24, 2, 20, 1000)

library(ggpubr)
ggarrange(p1, p2,
  ncol = 2, labels = "AUTO", common.legend = TRUE,
  font.label = list(size = 28, color = "black", face = "bold", family = NULL)
)
ggsave("plots/figS45.pdf", width = 16, height = 8, device = cairo_pdf)
