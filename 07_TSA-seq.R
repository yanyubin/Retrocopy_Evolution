##### TSA-seq and related data analysis #####

##### Plots #####

library(tidyverse)
library(ggsignif)
source("theme_publication.R")
source("format_pvalues.R")

factors <- c(
  "colocal_retrocopies", "noncolocal_retrocopies",
  "colocal_parents", "noncolocal_parents"
)


tsa_plot <- function(data, y, title, s) {
  ggplot(data, aes(x = Coloc, y = Score)) +
    geom_violin(aes(color = Coloc, fill = Coloc), scale = 1.6, alpha = 1.0) +
    geom_boxplot(width = 0.2, outlier.shape = NA, size = 1.5) +
    facet_grid(. ~ Type) +
    geom_signif(
      comparisons = list(
        c("Colocalized", "Noncolocalized")
      ),
      map_signif_level = TRUE, size = 0.8, textsize = 10, vjust = 0.5
    ) +
    ylim(min(data$Score), max(data$Score) * s) +
    labs(x = "", y = y, tag = "", title = title) +
    theme_publication(baseSize = 24) +
    theme(
      legend.position = "top",
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.title = element_blank()
    ) +
    scale_fill_discrete(labels = c("Colocalized", "Noncolocalized"))
}


read_tsa <- function(dataFile) {
  df <- read_tsv(dataFile, col_names = c("Type", "ID", "Score")) %>%
    separate_wider_delim(Type, "_", names = c("Coloc", "Type")) %>%
    mutate(Coloc = str_replace(Coloc, "colocal", "Colocalized")) %>%
    mutate(Coloc = str_replace(Coloc, "nonColocalized", "Noncolocalized")) %>%
    mutate(Type = str_replace(Type, "retrocopies", "Retrocopy")) %>%
    mutate(Type = str_replace(Type, "parents", "Parental gene"))

  df$Type <- factor(df$Type, levels = c("Retrocopy", "Parental gene"), ordered = TRUE)

  return(df)
}


##### SON TSA-seq #####
SON_TSA <- read_tsa("data/TSA-seq/K562_SON_TSA-Seq.txt")
tsa_plot(
  SON_TSA, bquote(bold(log[2] * "(pull-down/input)")),
  "K562 SON TSA-seq", 1.15
)

ggsave("plots/fig2C.pdf", width = 8, height = 8, device = cairo_pdf)


##### pSC35 TSA-Seq #####
pSC35_TSA <- read_tsa("data/TSA-seq/K562_pSC35_TSA-Seq.txt")
p1 <- tsa_plot(
  pSC35_TSA, bquote(bold(log[2] * "(pull-down/input)")),
  "K562 pSC35_TSA TSA-seq", 1.1
)

##### LaminAC TSA-Seq #####
laminAC_TSA <- read_tsa("data/TSA-seq/K562_LaminAC_TSA-Seq.txt")
p2 <- tsa_plot(
  laminAC_TSA, bquote(bold(log[2] * "(pull-down/input)")),
  "K562 lamin A/C TSA-seq", 1.15
)

##### LaminB TSA-Seq #####
laminB_TSA <- read_tsa("data/TSA-seq/K562_LaminB_TSA-Seq.txt")
p3 <- tsa_plot(
  laminB_TSA, bquote(bold(log[2] * "(pull-down/input)")),
  "K562 lamin B TSA-seq", 1.2
)

##### Speckle_distance #####
speckle_distance <- read_tsa("data/TSA-seq/K562_Speckle_distance.txt") %>%
  filter(Score < 1.35)

p4 <- tsa_plot(
  speckle_distance, bquote("Distance (μm)"), "K562 speckle distance",
  1.1
)

# Save plot
library(ggpubr)
ggarrange(p1, p2, p3, p4,
  labels = "AUTO", common.legend = TRUE,
  font.label = list(size = 28, color = "black", face = "bold", family = NULL)
)
ggsave("plots/figS8.pdf", width = 16, height = 18, device = cairo_pdf)


##### SON TSA Deciles #####
SON_Decile <- read_tsv("data/TSA-seq/K562_SON_Decile.txt",
  col_names = c("Type", "ID", "Group")
)

factors2 <- c(
  "Group_1", "Group_2", "Group_3", "Group_4", "Group_5",
  "Group_6", "Group_7", "Group_8", "Group_9", "Group_10"
)

SON_Decile$Group <- factor(SON_Decile$Group, levels = factors2, ordered = TRUE)

ggplot(SON_Decile, aes(x = Type, fill = Group)) +
  geom_bar(position = "fill") +
  labs(x = "", y = "Proportion", title = "K562") +
  theme_publication(baseSize = 28) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  theme(axis.title = element_text(size = 22)) +
  scale_x_discrete(labels = c(
    "Coloc_parent", "Coloc_retro",
    "Noncoloc_parent", "Noncoloc_retro"
  )) +
  scale_fill_discrete(
    name = "SON TSA\ndecile",
    labels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")
  ) +
  theme(legend.text.align = 0)

ggsave("plots/fig2D.pdf", width = 8, height = 8)



##### GRO-Seq #####
GRO_seq <- read_tsa("data/TSA-seq/K562_GRO-seq.txt")
GRO_seq$Score <- log2(GRO_seq$Score + 1)
p1 <- tsa_plot(
  GRO_seq, bquote(bold(log[2] * "(#reads + 1)")),
  "K562 GRO-seq", 1.1
)

##### Pol2 TSA-Seq #####
Pol2_TSA <- read_tsa("data/TSA-seq/K562_Pol2_TSA-Seq.txt")
p2 <- tsa_plot(
  Pol2_TSA, bquote(bold(log[2] * "(pull-down/input)")),
  "K562 Pol II TSA-seq", 1.15
)

# Save plot
library(ggpubr)
ggarrange(p1, p2,
  common.legend = TRUE,
  font.label = list(size = 28, color = "black", face = "bold", family = NULL)
)
ggsave("plots/fig3B.pdf", width = 14, height = 8)



##### Simulation that controls speckle distance #####
load("Rdata/3-simulations.Rdata")
library(tidyverse)
library(GenomicInteractions)

# Assign speckle distances to retrocopies and parental genes
dist1 <- read_tsv("data/TSA-seq/K562_Speckle_distance.txt",
  col_names = c("Type", "retroID", "Score")
)
human_retro <- human_retro %>% left_join(dist1, by = "retroID")
dist1 <- dist1 %>% dplyr::rename(parentID = retroID)
human_retro <- human_retro %>% left_join(dist1, by = "parentID")

# Read speckle distances for all protein coding genes
dist2 <- read_tsv("data/TSA-seq/human_coding_genes_Speckle_distance.txt",
  col_names = c("Type", "parentID", "Score")
) %>%
  left_join(human_pcgenes, by = "parentID")


humanCodingGenes <- GRanges(
  seqnames = dist2$parentChr,
  ranges = IRanges(
    start = dist2$parentStart,
    end = dist2$parentEnd
  ),
  strand = dist2$parentStrand,
  score = dist2$Score,
  names = dist2$parentID
)

humanRetro <- df2grange(human_retro, type = "retrocopy")


# Run simulation
library(doParallel)

# Setup backend to use many processors
totalCores <- detectCores()

# Leave one core to avoid overload your computer
cluster <- makeCluster(totalCores[1] - 1)
registerDoParallel(cluster)


x <- foreach(i = 1:1000, .combine = rbind) %dopar% {
  require(GenomicInteractions)
  require(tidyverse)
  require(plyranges)
  ParentTemp <- GRanges()
  for (j in 1:nrow(human_retro)) {
    d <- human_retro[j, ]$Score.y
    g <- plyranges::filter(humanCodingGenes, score > (d - 0.01) & score < (d + 0.01))
    Temp <- sample(g, 1)
    ParentTemp <- append(ParentTemp, Temp)
  }

  shuffleRP <- GenomicInteractions(humanRetro, ParentTemp, counts = 1)
  paired_hits <- findOverlaps(shuffleRP, K562_HiC)
  r <- length(unique(queryHits(paired_hits)))

  return(r)
}

# Stop cluster
stopCluster(cluster)

sim_df <- as.data.frame(x) %>%
  dplyr::rename(K562 = V1) %>%
  dplyr::mutate(Simulation = "Control speckle distance")

write_tsv(sim_df, "results/simulation_speckle_distance.tsv")


# Plot
library(tidyverse)
source("theme_publication.R")

df <- human_retro %>% dplyr::filter(Score.x < 1.5) # remove an outlier

p1 <- ggplot(df, aes(x = Score.x, y = Score.y, color = K562)) +
  geom_point(alpha = 0.3) +
  geom_density_2d() +
  labs(
    x = "Distance to speckle (μm)\nfor retrocopies ",
    y = "Distance to speckle (μm)\nfor parental genes "
  ) +
  theme_publication(baseSize = 24) +
  theme(legend.position = "top")

p2 <- sim_plot(human_retro, "K562", sim_df, 24, 2, 25, 1000) +
  theme(legend.position = "none")


# Save plot
library(ggpubr)
ggarrange(p1, p2,
  labels = "AUTO",
  font.label = list(size = 28, color = "black", face = "bold", family = NULL)
)
ggsave("plots/figS10.pdf", width = 16, height = 8, device = cairo_pdf)
