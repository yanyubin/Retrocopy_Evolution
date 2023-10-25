##### Function and conservation analysis in humans #####

#####  Ribo-seq RPKM #####
load("Rdata/5-expression.RData")

library(tidyverse)
library(ggsignif)
source("theme_publication.R")

human_retro$colocNum <- apply(human_retro, 1, function(x) length(which(x == "Colocalized")))
human_retro <- human_retro %>%
  dplyr::mutate(colocType = ifelse(colocNum >= 2, "≥2", as.character(colocNum)))


###
box_violin_plot <- function(df, yVar, size, y1, y2, y3) {
  require(ggplot2)
  require(ggsignif)

  ggplot(df, aes_string(x = "colocType", y = yVar)) +
    geom_violin(aes(color = colocType, fill = colocType), scale = 1.6, alpha = 0.5) +
    geom_boxplot(width = 0.2, outlier.shape = NA, size = size) +
    geom_signif(
      comparisons = list(
        c("0", "1"),
        c("1", "≥2"),
        c("0", "≥2")
      ), y_position = c(y1, y2, y3),
      map_signif_level = TRUE, size = 0.8, textsize = 10
    ) +
    labs(
      x = "Number of colocalized cell lines",
      # y = bquote(bold(log[10] * (RPKM))) # for Ribo-seq
      y = "LINSIGHT score" # for conservation score
    ) +
    theme_publication(baseSize = 24) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(size = 24)
    )
}

###
box_violin_plot2 <- function(df, yVar) {
  ggplot(df, aes_string(x = "colocal", y = yVar)) +
    geom_violin(aes(fill = factor(colocal))) +
    geom_boxplot(width = 0.3, outlier.shape = NA, size = 0.8) +
    facet_grid(. ~ cell_line) +
    geom_signif(
      comparisons = list(c("Colocalized", "Noncolocalized")),
      map_signif_level = TRUE, size = 0.8, textsize = 6
    ) +
    labs(x = "", y = "LINSIGHT score") +
    # labs(x = "", y = bquote(bold(log[10] * (RPKM)))) +
    guides(fill = guide_legend(title = NULL)) +
    theme_publication(baseSize = 24) +
    theme(
      legend.position = "top",
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    ) +
    scale_fill_discrete(labels = c("Colocalized", "Noncolocalized"))
}


###
retro <- read_tsv("data/RPFdb/human/human_retro_id.txt")

ribo_files <- list.files("data/RPFdb/human/RPKM/", pattern = NULL, full.names = TRUE)

for (f in ribo_files) {
  ribo <- read_tsv(f) %>%
    dplyr::select(-Gene_Name) %>%
    filter(!str_detect(Gene_ID, "_PAR_Y"))

  ribo$Gene_ID <- gsub("\\..*", "", ribo$Gene_ID)

  retro <- retro %>% left_join(ribo, by = "Gene_ID")
}

retro <- retro %>%
  rowwise() %>%
  mutate(
    RPKM_max = max(c_across(ERX432359:SRX2343261)),
    RPKM_mean = mean(c_across(ERX432359:SRX2343261), na.rm = TRUE)
  ) %>%
  dplyr::select(Gene_ID, annotatedRetroID, RPKM_max, RPKM_mean)

ribo_df <- human_retro %>%
  left_join(retro, by = "annotatedRetroID") %>%
  filter(RPKM_max >= 1)

ribo_df$colocType <- factor(ribo_df$colocType, levels = c("0", "1", "≥2"))

# Plot
box_violin_plot(ribo_df, "log10(RPKM_max)", 1.5, 3.6, 3.9, 4.2)
ggsave("plots/fig5A.pdf", width = 8, height = 10, device = cairo_pdf)


###
library(reshape2)

ribo_df2 <- human_retro %>% select(c(retroID, GM12878, HMEC, HUVEC, IMR90, K562, KBM7, NHEK))
ribo_df2 <- melt(ribo_df2,
  id.vars = "retroID",
  variable.name = "cell_line", value.name = "colocal"
)

ribo_df2 <- ribo_df2 %>%
  left_join(ribo_df, by = "retroID") %>%
  filter(!is.na(annotatedRetroID))

# Plot
box_violin_plot2(ribo_df2, "log10(RPKM_max)")
ggsave("plots/figS26.pdf", width = 12, height = 8)


##### Ribo-seq ORFs #####

orf <- read_tsv("data/RPFdb/human/human_ORFs.txt", col_names = F)
names(orf) <- c("orfID", "orfLength")
orf <- orf %>% mutate(transcriptID = str_sub(orfID, 1, 15))

n_orfs <- as.data.frame(table(orf$transcriptID)) %>% rename(annotatedRetroID = Var1)

uniq_orf <- unique(orf)

l_orfs <- uniq_orf %>%
  group_by(transcriptID) %>%
  summarise_at(vars(orfLength), list(meanLength = mean, maxLength = max)) %>%
  rename(annotatedRetroID = transcriptID)

orf_df <- human_retro %>%
  left_join(n_orfs, by = "annotatedRetroID") %>%
  left_join(l_orfs, by = "annotatedRetroID")

ribo_df <- ribo_df %>%
  left_join(orf_df[, c("retroID", "Freq")], by = "retroID") %>%
  filter(Freq > 0)

p1 <- box_violin_plot(ribo_df, "log10(RPKM_max)", 1.5, 3.6, 3.9, 4.2)

#
library(reshape2)
ribo_df2 <- human_retro %>% select(c(retroID, GM12878, HMEC, HUVEC, IMR90, K562, KBM7, NHEK))
ribo_df2 <- melt(ribo_df2,
  id.vars = "retroID",
  variable.name = "cell_line", value.name = "colocal"
)
ribo_df2 <- ribo_df2 %>%
  left_join(ribo_df, by = "retroID") %>%
  filter(!is.na(annotatedRetroID))

p2 <- box_violin_plot2(ribo_df2, "log10(RPKM_max)")

# Save plot
library(ggpubr)
ggarrange(p1, p2,
  widths = c(1, 2),
  ncol = 2, labels = "AUTO",
  font.label = list(size = 28, color = "black", face = "bold", family = NULL)
)
ggsave("plots/figS27.pdf", width = 20, height = 10, device = cairo_pdf)



##### PhastCons score #####

phast <- read_tsv("data/Conservation_Score/human_retro.PhastCons46way.tab")
df <- human_retro %>%
  left_join(phast, by = "retroID")
df$colocType <- factor(df$colocType, levels = c("0", "1", "≥2"))

# Save plot
box_violin_plot(df, "mean", 1.5, 1.0, 1.08, 1.16)
ggsave("plots/fig5B.pdf", width = 8, height = 10, device = cairo_pdf)


###
library(reshape2)
df <- human_retro %>% dplyr::select(c(retroID, GM12878, HMEC, HUVEC, IMR90, K562, KBM7, NHEK))
df <- melt(df, id.vars = "retroID", variable.name = "cell_line", value.name = "colocal")
df <- df %>% left_join(phast, by = "retroID")

# Save plot
box_violin_plot2(df, "mean")
ggsave("plots/figS28.pdf", width = 12, height = 8)


##### LINSIGHT score #####

linsight <- read_tsv("data/Conservation_Score/human_retro.linsight.tab")
df <- human_retro %>%
  left_join(linsight, by = "retroID")
df$colocType <- factor(df$colocType, levels = c("0", "1", "≥2"))

p1 <- box_violin_plot(df, "mean", 1.5, 0.95, 1.03, 1.1)


###
df <- human_retro %>% dplyr::select(c(retroID, GM12878, HMEC, HUVEC, IMR90, K562, KBM7, NHEK))
df <- melt(df, id.vars = "retroID", variable.name = "cell_line", value.name = "colocal")
df <- df %>% left_join(linsight, by = "retroID")

p2 <- box_violin_plot2(df, "mean")

library(ggpubr)
ggarrange(p1, p2,
  widths = c(1, 2),
  ncol = 2, labels = "AUTO",
  font.label = list(size = 28, color = "black", face = "bold", family = NULL)
)
ggsave("plots/figS29.pdf", width = 20, height = 10, device = cairo_pdf)



##### GWAS hits #####

library(tidyverse)
library(ggsignif)
source("theme_publication.R")

load("Rdata/2-colocalization_250kb.Rdata")

gwas <- read_table("data/Conservation_Score/human_retro_inter_gwas_promoter1k.bed")

human_retro <- human_retro %>% left_join(gwas, by = "retroID")

df1 <- human_retro %>%
  mutate(GWAS_hits = ifelse(num_overlap == 0, "No", "Yes")) %>%
  group_by(GM12878, GWAS_hits) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  rename(coloc = GM12878) %>%
  cbind(data.frame(cell_line = rep("GM12878", 4)))

df2 <- human_retro %>%
  mutate(GWAS_hits = ifelse(num_overlap == 0, "No", "Yes")) %>%
  group_by(K562, GWAS_hits) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  rename(coloc = K562) %>%
  cbind(data.frame(cell_line = rep("K562", 4)))

cell_lines <- c("GM12878", "HMEC", "HUVEC", "IMR90", "K562", "KBM7", "NHEK")

for (cell_line in cell_lines) {
  df <- human_retro %>%
    mutate(GWAS_hits = ifelse(num_overlap == 0, "No", "Yes")) %>%
    group_by(!!!rlang::syms(cell_line), GWAS_hits) %>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n)) %>%
    rename(coloc = cell_line) %>%
    cbind(data.frame(cell_line = rep(cell_line, 4)))

  assign(cell_line, df)
}

df <- rbind(GM12878, HMEC, HUVEC, IMR90, K562, KBM7, NHEK)
df$GWAS_hits <- factor(df$GWAS_hits, levels = c("Yes", "No"))

fisher.test(matrix(GM12878$n, ncol = 2)) # p-value = 0.0002124
fisher.test(matrix(HMEC$n, ncol = 2)) # p-value = 2.977e-05
fisher.test(matrix(HUVEC$n, ncol = 2)) # p-value = 2.099e-05
fisher.test(matrix(IMR90$n, ncol = 2)) # p-value = 0.00916
fisher.test(matrix(K562$n, ncol = 2)) # p-value = 8.406e-05
fisher.test(matrix(KBM7$n, ncol = 2)) # p-value = 0.007203
fisher.test(matrix(NHEK$n, ncol = 2)) # p-value = 0.0603


ggplot(df, aes(x = coloc, y = freq, fill = GWAS_hits)) +
  geom_bar(stat = "identity") +
  facet_grid(. ~ cell_line) +
  geom_signif(
    comparisons = list(c("Colocalized", "Noncolocalized")), y_position = c(1.05),
    tip_length = c(0, 0),
    annotations = c("***"), size = 1.1, textsize = 8
  ) +
  geom_text(aes(label = paste0(sprintf("%1.1f", freq * 100), "%", "\n(", n, ")")),
    position = position_stack(vjust = 0.5), size = 6
  ) +
  ylim(c(0, 1.2)) +
  labs(x = "", y = "Percentage") +
  theme_publication(baseSize = 28) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  theme(axis.title = element_text(size = 22)) +
  scale_x_discrete(labels = c("Colocalized", "Noncolocalized")) +
  theme(legend.position = "top")

ggsave("plots/fig5C.pdf", width = 18, height = 12)



##### Correlation between PhastConse score and expression/translation #####

df <- ribo_df %>% left_join(phast, by = "retroID")

exp <- df %>% dplyr::select(ends_with("Pap.R"))
exp$exp_max <- apply(exp, 1, max, na.rm = TRUE)
exp$exp_mean <- apply(exp, 1, mean, na.rm = TRUE)

df <- cbind(df, exp[c("exp_max", "exp_mean")])


library(reshape2)
df_m <- df %>%
  dplyr::select(exp_max, RPKM_max, mean) %>%
  melt(id.vars = "mean", variable.name = "type", value.name = "FPKM") %>%
  mutate(logFPKM = log2(FPKM + 1))

df_m$type <- gsub("exp_max", "RNA-seq", df_m$type)
df_m$type <- gsub("RPKM_max", "Ribo-seq", df_m$type)
df_m$type <- factor(df_m$type, levels = c("RNA-seq", "Ribo-seq"), ordered = TRUE)

library("ggpubr")
ggscatter(df_m,
  x = "mean", y = "logFPKM", facet.by = "type", size = 1,
  add = "reg.line", add.params = list(size = 1.5),
  conf.int = TRUE, color = "type", palette = "npg",
  cor.coef = TRUE, cor.method = "pearson",
  xlab = "PhastCons score", ylab = "log2(RPKM + 1)",
  cor.coef.size = 6, alpha = 0.3
) +
  theme(
    text = element_text(size = 16),
    plot.title = element_text(size = 16, face = "bold", colour = "black"),
    axis.title.x = element_text(size = 16, face = "bold", colour = "black"),
    axis.title.y = element_text(size = 16, face = "bold", colour = "black"),
  ) +
  theme(legend.position = "none")

ggsave("plots/fig5D.pdf", width = 9, height = 6)
