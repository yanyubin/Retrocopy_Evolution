##### Mouse expression and conservation data analysis #####

library(tidyverse)
load("Rdata/2-colocalization_250kb.Rdata")
library(ggpubr)
source("theme_publication.R")


##### RNA-seq #####
rep1.R <- read_tsv("data/ENCODE_RNA_Dashboard/CH12.LX/rep1_mm10.tsv") %>%
  rename(annotatedRetroID = `Gene ID`)
rep1.R$annotatedRetroID <- gsub("\\..*", "", rep1.R$annotatedRetroID)

rep2.R <- read_tsv("data/ENCODE_RNA_Dashboard/CH12.LX/rep2_mm10.tsv") %>%
  rename(annotatedRetroID = `Gene ID`)
rep2.R$annotatedRetroID <- gsub("\\..*", "", rep2.R$annotatedRetroID)


rep1.P <- read_tsv("data/ENCODE_RNA_Dashboard/CH12.LX/rep1_mm10.tsv") %>%
  rename(parentID = `Gene ID`)
rep1.P$parentID <- gsub("\\..*", "", rep1.P$parentID)

rep2.P <- read_tsv("data/ENCODE_RNA_Dashboard/CH12.LX/rep2_mm10.tsv") %>%
  rename(parentID = `Gene ID`)
rep2.P$parentID <- gsub("\\..*", "", rep2.P$parentID)


df <- mouse_retro %>%
  left_join(rep1.R, by = "annotatedRetroID") %>%
  left_join(rep2.R, by = "annotatedRetroID") %>%
  left_join(rep1.P, by = "parentID") %>%
  left_join(rep2.P, by = "parentID") %>%
  rowwise() %>%
  mutate(FPKM.R = mean(c(FPKM.x, FPKM.y))) %>%
  mutate(FPKM.P = mean(c(FPKM.x.x, FPKM.y.y))) %>%
  dplyr::select(retroID:CH12.LX, FPKM.R, FPKM.P) %>%
  distinct(retroID, .keep_all = TRUE)


C <- df %>% filter(CH12.LX == "Colocalized")
N <- df %>% filter(CH12.LX == "Noncolocalized")

t.R <- wilcox.test(log10(C$FPKM.R + 1), log10(N$FPKM.R + 1))
t.P <- wilcox.test(log10(C$FPKM.P + 1), log10(N$FPKM.P + 1))


p1 <- ggplot(df, aes(x = log10(FPKM.R + 1))) +
  stat_ecdf(aes_string(color = "CH12.LX"),
    geom = "step", size = 1.5
  ) +
  scale_color_manual(
    values = c(rgb(231 / 255, 126 / 255, 114 / 255), rgb(85 / 255, 188 / 255, 194 / 255)),
    name = "", labels = c("Colocalized", "Noncolocalized")
  ) +
  annotate("text", x = 1, y = 0.75, label = sprintf("p = %.1e", t.R$p.value), size = 8) +
  labs(
    x = bquote(bold(log[10] * "(FPKM + 1)")), y = "Cumulative fraction",
    title = "RNA-seq (Retrocopies)"
  ) +
  theme_publication(baseSize = 24)


p2 <- ggplot(df, aes(x = log10(FPKM.P + 1))) +
  stat_ecdf(aes_string(color = "CH12.LX"),
    geom = "step", size = 1.5
  ) +
  scale_color_manual(
    values = c(rgb(231 / 255, 126 / 255, 114 / 255), rgb(85 / 255, 188 / 255, 194 / 255)),
    name = "", labels = c("Colocalized", "Noncolocalized")
  ) +
  annotate("text", x = 1.18, y = 1, label = sprintf("p = %.1e", t.P$p.value), size = 8) +
  labs(
    x = bquote(bold(log[10] * "(FPKM + 1)")), y = "Cumulative fraction",
    title = "RNA-seq (Parental genes)"
  ) +
  theme_publication(baseSize = 24)

library(ggpubr)
ggarrange(p1, p2, common.legend = TRUE)
ggsave("plots/figS34.pdf", width = 12, height = 6)



##### Ribo-seq #####

retro <- read_tsv("data/RPFdb/mouse/mouse_retro_id.txt")

Files <- list.files("data/RPFdb/mouse/RPKM/", pattern = NULL, full.names = TRUE)

for (f in Files) {
  ribo <- read_tsv(f) %>%
    dplyr::select(-Gene_Name) %>%
    filter(!str_detect(Gene_ID, "_PAR_Y"))

  ribo$Gene_ID <- gsub("\\..*", "", ribo$Gene_ID)

  retro <- retro %>% left_join(ribo, by = "Gene_ID")
}


retro <- retro %>%
  rowwise() %>%
  mutate(
    RPKM_max = max(c_across(SRX2943749:SRX172394)),
    RPKM_mean = mean(c_across(SRX2943749:SRX172394), na.rm = TRUE)
  ) %>%
  dplyr::select(Gene_ID, annotatedRetroID, RPKM_max, RPKM_mean)


df <- df %>% left_join(retro, by = "annotatedRetroID")

ribo_df <- df %>%
  filter(RPKM_max >= 1)


# plot
library(ggsignif)

p3 <- ggplot(ribo_df, aes(x = `CH12.LX`, y = log10(RPKM_max))) +
  geom_violin(aes(fill = factor(`CH12.LX`))) +
  geom_boxplot(width = 0.3, outlier.shape = NA, size = 1) +
  geom_signif(
    comparisons = list(c("Colocalized", "Noncolocalized")),
    map_signif_level = TRUE, size = 0.8, textsize = 8
  ) +
  labs(x = "", y = bquote(bold(log[10] * (RPKM))), title = "Ribo-seq") +
  guides(fill = guide_legend(title = NULL)) +
  theme_publication(baseSize = 24) +
  theme(
    legend.position = "top",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  scale_fill_discrete(labels = c("Colocalized", "Noncolocalized"))



##### PhastCons score #####

phast <- read_tsv("data/Conservation_Score/mm10.60way.phastCons.tab")

df <- df %>%
  left_join(phast, by = "retroID")

p4 <- ggplot(df, aes(x = `CH12.LX`, y = mean)) +
  geom_violin(aes(fill = factor(`CH12.LX`))) +
  geom_boxplot(width = 0.15, outlier.shape = NA, size = 1) +
  geom_signif(
    comparisons = list(c("Colocalized", "Noncolocalized")),
    map_signif_level = TRUE, size = 0.8, textsize = 8
  ) +
  labs(x = "", y = "PhastCons score", title = "Conservation") +
  guides(fill = guide_legend(title = NULL)) +
  theme_publication(baseSize = 24) +
  theme(
    legend.position = "top",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  scale_fill_discrete(labels = c("Colocalized", "Noncolocalized"))


ggarrange(p3, p4, common.legend = TRUE)
ggsave("plots/figS35.pdf", width = 12, height = 10)



##### Correlation #####

library(reshape2)
df_m <- df %>%
  dplyr::select(FPKM.R, RPKM_max, mean) %>%
  melt(id.vars = "mean", variable.name = "type", value.name = "FPKM")

df_m <- df_m %>% mutate(logFPKM = log2(FPKM + 1))
df_m$type <- gsub("FPKM.R", "RNA-seq", df_m$type)
df_m$type <- gsub("RPKM_max", "Ribo-seq", df_m$type)

df_m$type <- factor(df_m$type, levels = c("RNA-seq", "Ribo-seq"), ordered = TRUE)

library("ggpubr")
ggscatter(df_m,
  x = "mean", y = "logFPKM", facet.by = "type", size = 1,
  add = "reg.line", add.params = list(size = 1.5),
  conf.int = TRUE, color = "type", palette = "npg",
  cor.coef = TRUE, cor.method = "pearson",
  xlab = "PhastCons score", ylab = "log2(RPKM+1)",
  cor.coef.size = 6, alpha = 0.3
) +
  theme(
    text = element_text(size = 16),
    plot.title = element_text(size = 16, face = "bold", colour = "black"),
    axis.title.x = element_text(size = 16, face = "bold", colour = "black"),
    axis.title.y = element_text(size = 16, face = "bold", colour = "black"),
  ) +
  theme(legend.position = "none")

ggsave("./figures/figS36.pdf", width = 8, height = 6)
