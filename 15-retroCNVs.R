##### retroCNV analysis #####

library(tidyverse)
source("theme_publication.R")

retroCNV <- read_table("data/retroCNVs/tableS3.tsv") %>%
  select(retroChr:allele_count, -parentID)

##### SNP derived allele frequencies #####

snps <- read_tsv("data/retroCNVs/5pop_SNPs_ac.frq.count", col_names = FALSE, skip = 1)
names(snps) <- c("chr", "pos", "n_allele", "n_ind", "ac1", "allele_count")

snps <- snps %>%
  filter(allele_count != 0) %>%
  mutate(AC = ifelse(allele_count > 100, 51, allele_count / 2))

breaks <- seq(0, 55, 5)

snps$AC2 <- cut(snps$AC, breaks = breaks)

snp_cut <- snps %>%
  group_by(AC2) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))


##### RetroCNV insertion frequencies #####
retroCNV$colocNum <- apply(
  X = retroCNV, MARGIN = 1,
  function(t) {
    sum(grepl(pattern = "C", x = t, fixed = TRUE))
  }
)

retroCNV <- retroCNV %>%
  mutate(colocType = ifelse(colocNum >= 2, "≥2", as.character(colocNum)))

retroCNV <- retroCNV %>%
  mutate(AC = ifelse(allele_count > 50, 51, allele_count))

breaks <- seq(0, 55, 5)

retroCNV$AC2 <- cut(retroCNV$AC, breaks = breaks)


retroCNV_cut1 <- retroCNV %>%
  group_by(AC2) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))


snp_cut$type <- c(rep("SNPs", 11))
retroCNV_cut1$type <- c(rep("retroCNVs", 11))

df <- rbind(snp_cut, retroCNV_cut1)

df <- mutate(df, AC2 = fct_recode(AC2, "(50,491]" = "(50,55]"))

ggplot(data = df, aes(x = AC2, y = freq, fill = type)) +
  geom_bar(stat = "identity", position = "dodge") +
  ylab("Frequency of retroCNVs") +
  xlab("Number of individuals") +
  theme_publication(baseSize = 20) +
  guides(fill = guide_legend(title = NULL)) +
  theme(legend.position = "top", axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("plots/figS38.pdf", width = 9.7, height = 6)


##### SFS for retrocopies colocalized in 0, 1, and >= 2 cell lines #####
retroCNV_cut2 <- retroCNV %>%
  group_by(colocType, AC2) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))

retroCNV_cut2 <- mutate(retroCNV_cut2, AC2 = fct_recode(AC2, "(50,491]" = "(50,55]"))
retroCNV_cut2$colocType <- factor(retroCNV_cut2$colocType, levels = c("0", "1", "≥2"))

ggplot(data = retroCNV_cut2, aes(x = AC2, y = freq, fill = colocType)) +
  geom_bar(stat = "identity", position = "dodge") +
  ylab("Frequency of retroCNVs") +
  xlab("Number of individuals") +
  theme_publication(baseSize = 20) +
  guides(fill = guide_legend(title = "Number of colocalized cell lines")) +
  theme(legend.position = "top", axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("plots/fig6B.pdf", width = 9.7, height = 6, device = cairo_pdf)


##### SFS for retrocopies colocalized in individual cell lines #####
retroCNV_cut2 <- retroCNV %>%
  group_by(NHEK, AC2) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))

retroCNV_cut2 <- mutate(retroCNV_cut2, AC2 = fct_recode(AC2, "(50,491]" = "(50,55]"))

p7 <- ggplot(data = retroCNV_cut2, aes(x = AC2, y = freq, fill = NHEK)) +
  geom_bar(stat = "identity", position = "dodge") +
  ylab("Frequency of retroCNVs") +
  xlab("Number of individuals") +
  theme_publication(baseSize = 20) +
  labs(title = "NHEK") +
  guides(fill = guide_legend(title = NULL)) +
  theme(legend.position = "top", axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_discrete(labels = c("Colocalized", "Noncolocalized"))


library(ggpubr)
ggarrange(p1, p2, p3, p4, p5, p6, p7,
  ncol = 4, nrow = 2,
  common.legend = TRUE,
  font.label = list(color = "black", size = 18)
)

ggsave("plots/figS41.pdf", width = 20, height = 10)




##### Density of conserved sites #####
cs <- read_tsv("data/retroCNVs/retroCNV_hg19_100kb_cs_density.bed",
  col_names = FALSE
) %>% rename(cs_density = X8)

retroCNV <- cbind(retroCNV, cs[, "cs_density"])
retroCNV$colocType <- factor(retroCNV$colocType, levels = c("0", "1", "≥2"))

library(ggsignif)
p1 <- ggplot(retroCNV, aes(x = colocType, y = cs_density)) +
  geom_violin(aes(color = colocType, fill = colocType), scale = 1.6, alpha = 0.5) +
  geom_boxplot(width = 0.2, outlier.shape = NA, size = 1.5) +
  geom_signif(
    comparisons = list(
      c("0", "1"),
      c("1", "≥2"),
      c("0", "≥2")
    ), y_position = c(0.32, 0.36, 0.4),
    map_signif_level = TRUE, size = 0.8, textsize = 8
  ) +
  labs(
    x = "Number of colocalized cell lines", y = "Density of conserved sites",
    tag = ""
  ) +
  theme_publication(baseSize = 24) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 24)
  )


###
library(reshape2)
df <- retroCNV %>% dplyr::select(c(retroEnd, GM12878, HMEC, HUVEC, IMR90, K562, KBM7, NHEK))
df <- melt(df, id.vars = "retroEnd", variable.name = "cell_line", value.name = "colocal")

df <- df %>% left_join(retroCNV, by = "retroEnd")

p2 <- ggplot(df, aes(x = colocal, y = cs_density)) +
  geom_violin(aes(fill = factor(colocal))) +
  geom_boxplot(width = 0.3, outlier.shape = NA, size = 0.8) +
  facet_grid(. ~ cell_line) +
  geom_signif(
    comparisons = list(c("C", "N")),
    map_signif_level = TRUE, size = 0.8, textsize = 8
  ) +
  labs(x = "", y = "Density of conserved sites", tag = "") +
  guides(fill = guide_legend(title = NULL)) +
  theme_publication(baseSize = 24) +
  theme(
    legend.position = "top",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  scale_fill_discrete(labels = c("Colocalized", "Noncolocalized"))


library(ggpubr)
ggarrange(p1, p2,
  widths = c(1, 2),
  ncol = 2, labels = "AUTO",
  font.label = list(size = 28, color = "black", face = "bold", family = NULL)
)
ggsave("plots/figS42.pdf", width = 20, height = 10, device = cairo_pdf)



##### DAF of flanking sites #####

###
breaks <- seq(0, 105, 5)

###
f0 <- "data/retroCNVs/retroCNV_hg19_C0_25000.frq.count"
snps0 <- read_tsv(f0, col_names = FALSE, skip = 1)
names(snps0) <- c("chr", "pos", "n_allele", "n_ind", "ac1", "allele_count")

snps0 <- snps0 %>%
  filter(allele_count != 0) %>%
  mutate(AC = ifelse(allele_count > 100, 101, allele_count))

snps0$AC2 <- cut(snps0$AC, breaks = breaks)

snp_cut0 <- snps0 %>%
  group_by(AC2) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))

snp_cut0$type <- c(rep("0", 21))


###
f1 <- "data/retroCNVs/retroCNV_hg19_C1_25000.frq.count"
snps1 <- read_tsv(f1, col_names = FALSE, skip = 1)
names(snps1) <- c("chr", "pos", "n_allele", "n_ind", "ac1", "allele_count")

snps1 <- snps1 %>%
  filter(allele_count != 0) %>%
  mutate(AC = ifelse(allele_count > 100, 101, allele_count))

snps1$AC2 <- cut(snps1$AC, breaks = breaks)

snp_cut1 <- snps1 %>%
  group_by(AC2) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))

snp_cut1$type <- c(rep("1", 21))


###
f2 <- "data/retroCNVs/retroCNV_hg19_C2_25000.frq.count"
snps2 <- read_tsv(f2, col_names = FALSE, skip = 1)
names(snps2) <- c("chr", "pos", "n_allele", "n_ind", "ac1", "allele_count")

snps2 <- snps2 %>%
  filter(allele_count != 0) %>%
  mutate(AC = ifelse(allele_count > 100, 101, allele_count))

snps2$AC2 <- cut(snps2$AC, breaks = breaks)

snp_cut2 <- snps2 %>%
  group_by(AC2) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))

snp_cut2$type <- c(rep("≥2", 21))


df <- rbind(snp_cut0, snp_cut1, snp_cut2)

df <- mutate(df, AC2 = fct_recode(AC2, "(100,981]" = "(100,105]"))
df$type <- factor(df$type, levels = c("0", "1", "≥2"))

ggplot(data = df, aes(x = AC2, y = freq, fill = type)) +
  geom_bar(stat = "identity", position = "dodge") +
  ylab("Frequency of SNPs") +
  xlab("Derived allele frequency") +
  theme_publication(baseSize = 20) +
  guides(fill = guide_legend(title = "Number of colocalized cell lines")) +
  theme(legend.position = "top", axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("plots/figS43.pdf", width = 9.7, height = 6, device = cairo_pdf)



########## Simulations for retroCNVs ##########

load("Rdata/3-simulations.Rdata")

# read data
library(tidyverse)
retroCNV_LF <- read_table("data/retroCNVs/tableS3.tsv") %>%
  dplyr::select(retroChr:allele_count) %>%
  dplyr::filter(allele_count == 1) %>%
  mutate(across(GM12878:NHEK, str_replace, "C", "Colocalized")) %>%
  mutate(across(GM12878:NHEK, str_replace, "N", "Noncolocalized"))


retroCNV_HF <- read_table("data/retroCNVs/tableS3.tsv") %>%
  dplyr::select(retroChr:allele_count) %>%
  dplyr::filter(allele_count >= 28) %>%
  mutate(across(GM12878:NHEK, str_replace, "C", "Colocalized")) %>%
  mutate(across(GM12878:NHEK, str_replace, "N", "Noncolocalized"))


# run simulations

library(doParallel)
# Setup backend to use many processors
totalCores <- detectCores()

# Leave one core to avoid overload your computer
cluster <- makeCluster(totalCores[1] - 4)
registerDoParallel(cluster)


n <- 1000

# Singleton retroCNVs
for (s in hic_data[1:7]) {
  assign(paste0("sim_", s), run_sim(
    get(paste0(s, "_HiC")), retroCNV_LF, hg19_chr_sizes,
    human_genes, human_pcgenes, n
  ))
}

# Save simulation results
sim <- cbind(sim_GM12878, sim_HMEC, sim_HUVEC, sim_IMR90, sim_K562, sim_KBM7, sim_NHEK)
names(sim) <- hic_data[1:7]

sim[, "Simulation"] <- as.vector(sapply(1:5, function(i) rep(i, n)))
readr::write_tsv(sim, "results/simulation_retroCNV_LF.tsv")


# High-frequency retroCNVs
for (s in hic_data[1:7]) {
  assign(paste0("sim_", s), run_sim(
    get(paste0(s, "_HiC")), retroCNV_HF, hg19_chr_sizes,
    human_genes, human_pcgenes, n
  ))
}

# Save simulation results
sim <- cbind(sim_GM12878, sim_HMEC, sim_HUVEC, sim_IMR90, sim_K562, sim_KBM7, sim_NHEK)
names(sim) <- hic_data[1:7]

sim[, "Simulation"] <- as.vector(sapply(1:5, function(i) rep(i, n)))
readr::write_tsv(sim, "results/simulation_retroCNV_HF.tsv")

# Plot LF
sim_LF <- read_tsv("results/simulation_retroCNV_LF.tsv")
sim_LF$Simulation <- as.factor(sim_LF$Simulation)

p1 <- sim_plot(retroCNV_LF, "GM12878", sim_LF, 24, 1.8, 4.2, 1000) +
  labs(title = "Singleton retroCNVs")
p2 <- sim_plot(retroCNV_LF, "HMEC", sim_LF, 16, 1, 5, 1000)
p3 <- sim_plot(retroCNV_LF, "HUVEC", sim_LF, 16, 1, 5, 1000)
p4 <- sim_plot(retroCNV_LF, "IMR90", sim_LF, 16, 1, 5, 1000)
p5 <- sim_plot(retroCNV_LF, "K562", sim_LF, 16, 1, 5, 1000)
p6 <- sim_plot(retroCNV_LF, "KBM7", sim_LF, 16, 1, 5, 1000)
p7 <- sim_plot(retroCNV_LF, "NHEK", sim_LF, 16, 1, 5, 1000)

# Save plot
library(ggpubr)
ggarrange(p2, p3, p4, p5, p6, p7,
  common.legend = TRUE, labels = "AUTO",
  font.label = list(color = "black", size = 18)
)

ggsave("plots/figS39.pdf", width = 12, height = 8)


# Plot HF
sim_HF <- read_tsv("results/simulation_retroCNV_HF.tsv")
sim_HF$Simulation <- as.factor(sim_HF$Simulation)

p11 <- sim_plot(retroCNV_HF, "GM12878", sim_HF, 24, 1.8, 5, 1000) +
  labs(title = "High-frequency retroCNVs")
p22 <- sim_plot(retroCNV_HF, "HMEC", sim_HF, 16, 1, 5, 1000)
p33 <- sim_plot(retroCNV_HF, "HUVEC", sim_HF, 16, 1, 5, 1000)
p44 <- sim_plot(retroCNV_HF, "IMR90", sim_HF, 16, 1, 5, 1000)
p55 <- sim_plot(retroCNV_HF, "K562", sim_HF, 16, 1, 5, 1000)
p66 <- sim_plot(retroCNV_HF, "KBM7", sim_HF, 16, 1, 5, 1000)
p77 <- sim_plot(retroCNV_HF, "NHEK", sim_HF, 16, 1, 5, 1000)

# Save plot HF
library(ggpubr)
ggarrange(p22, p33, p44, p55, p66, p77,
  common.legend = TRUE, labels = "AUTO",
  font.label = list(color = "black", size = 18)
)

ggsave("plots/figS40.pdf", width = 12, height = 8)


# Save plot GM12878
ggarrange(p1, p11, common.legend = TRUE)

ggsave("plots/fig6A.pdf", width = 12, height = 6)
