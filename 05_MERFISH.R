##### Multiplex FISH data analysis #####
load("Rdata/4-SPRITE.Rdata")

##### MERFISH (Su et al. 2020; IMR90) #####
library(reticulate)

# The pkl file is obtained using data and script of Su et al. (2020)
my_data <- py_load_object("data/interchrom_contacts/FISH/IMR90_MERFISH_contact_mat.pkl")

# Get distance and proximity frequency matrices (1041 loci)
df_dist <- as.data.frame(my_data[[1]])

df_prox_1um <- as.data.frame(my_data[[2]])[, 3124:4164]
names(df_prox_1um) <- c(paste0("V", 1:1041, sep = ""))


# Convert matrices to bedpe files
mat2bedpe <- function(mat) {
  require(tidyverse)

  sparse_mat <- full2sparse(as.matrix(mat))

  # read loci coordinates
  loci <- read_tsv("data/interchrom_contacts/FISH/IMR90_MERFISH_1041_loci_coords.txt",
    col_names = F
  ) %>%
    mutate(X2 = X2 - 700000, X3 = X3 + 700000) # expand to 1.5 Mb

  names(loci) <- c("chrom", "start", "end", "region1")
  mat_bedpe <- sparse_mat %>% left_join(loci, by = "region1")

  names(loci) <- c("chrom", "start", "end", "region2")
  mat_bedpe <- mat_bedpe %>% left_join(loci, by = "region2")

  mat_bedpe <- mat_bedpe %>%
    mutate(name = ".", strand1 = ".", strand2 = ".") %>%
    dplyr::select(chrom.x:end.y, "name", "IF", "strand1", "strand2") %>%
    dplyr::rename(distance = IF, `#chrom.x` = chrom.x)

  return(mat_bedpe)
}


dist_bedpe <- mat2bedpe(df_dist)
write_tsv(
  dist_bedpe,
  "data/interchrom_contacts/FISH/IMR90_MERFISH_1.5Mb_distance_mat.bedpe"
)

prox_bedpe <- mat2bedpe(df_prox_1um)
write_tsv(
  prox_bedpe,
  "data/interchrom_contacts/FISH/IMR90_MERFISH_1.5Mb_proximity_freq_mat.bedpe"
)



##### Comapare spatial distance #####
library(tidyverse)
human_retro_hg38 <- read_tsv("data/interchrom_contacts/FISH/human_retro_inter_hg38.bedpe") %>%
  dplyr::rename(retroChr = `#retroChr`)

library(GenomicInteractions)
rp_pairs_hg38 <- makeGenomicInteractionsFromFile(
  "data/interchrom_contacts/FISH/human_retro_inter_hg38.bedpe",
  type = "bedpe",
  experiment_name = "retro-parent pairs",
  description = "human"
)

merfish_dist <- makeGenomicInteractionsFromFile(
  "data/interchrom_contacts/FISH/IMR90_MERFISH_1.5Mb_distance_mat.bedpe",
  type = "bedpe",
  experiment_name = "merfish 1.5Mb",
  description = "IMR90"
)

merfish_prox <- makeGenomicInteractionsFromFile(
  "data/interchrom_contacts/FISH/IMR90_MERFISH_1.5Mb_proximity_freq_mat.bedpe",
  type = "bedpe",
  experiment_name = "merfish 1.5Mb",
  description = "IMR90"
)


m_dist <- as.data.frame(mergeByOverlaps(rp_pairs_hg38, merfish_dist, ignore.strand = TRUE)) %>%
  dplyr::select(rp_pairs_hg38.name, merfish_dist.counts) %>%
  dplyr::rename(retroID = rp_pairs_hg38.name, dist = merfish_dist.counts) %>%
  dplyr::group_by(retroID) %>%
  dplyr::summarise(dist_um = mean(dist))

m_prox <- as.data.frame(mergeByOverlaps(rp_pairs_hg38, merfish_prox, ignore.strand = TRUE)) %>%
  dplyr::select(rp_pairs_hg38.name, merfish_prox.counts) %>%
  dplyr::rename(retroID = rp_pairs_hg38.name, prox = merfish_prox.counts) %>%
  dplyr::group_by(retroID) %>%
  dplyr::summarise(prox_freq = mean(prox))


human_retro_hg38 <- human_retro_hg38 %>%
  left_join(m_dist, by = "retroID") %>%
  left_join(m_prox, by = "retroID")


# Plot
library(ggsignif)
source("format_pvalues.R")

# p = 3.5x10-6
p1 <- ggplot(human_retro_hg38, aes(x = IMR90, y = dist_um)) +
  geom_violin(aes(color = IMR90, fill = IMR90), scale = 1.6, alpha = 1.0) +
  geom_boxplot(width = 0.2, outlier.shape = NA, size = 1.5) +
  geom_signif(
    comparisons = list(
      c("Colocalized", "Noncolocalized")
    ),
    map_signif_level = TRUE, size = 1.1, textsize = 8
  ) +
  labs(x = "", y = "Spatial distance (μm)") +
  theme_publication(baseSize = 24) +
  theme(
    legend.position = "none",
    legend.title = element_blank()
  )

# p = 2.0x10-13
p2 <- ggplot(human_retro_hg38, aes(x = IMR90, y = prox_freq)) +
  geom_violin(aes(color = IMR90, fill = IMR90), scale = 1.6, alpha = 1.0) +
  geom_boxplot(width = 0.2, outlier.shape = NA, size = 1.5) +
  geom_signif(
    comparisons = list(
      c("Colocalized", "Noncolocalized")
    ),
    map_signif_level = TRUE, size = 1.1, textsize = 8
  ) +
  labs(x = "", y = "Proximity frequency") +
  theme_publication(baseSize = 24) +
  theme(
    legend.position = "none",
    legend.title = element_blank()
  )


library(ggpubr)
ggarrange(p1, p2, labels = "AUTO", font.label = list(size = 24, face = "bold"))

ggsave("plots/figS5.pdf", width = 12, height = 8, device = cairo_pdf)
