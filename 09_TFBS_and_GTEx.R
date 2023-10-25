##### TFBS data analysis #####

library(tidyverse)
library(ggsignif)
library(ggpubr)
library(reshape2)

source("theme_publication.R")
load("Rdata/5-expression.RData")


split_tibble <- function(tibble, column = "col") {
  tibble %>%
    split(., .[, column]) %>%
    lapply(., function(x) x[, setdiff(names(x), column)])
}


shared_tfbs_plot <- function(fimoFile) {
  require(vecsets)

  tfbs <- read_tsv(fimoFile)
  s <- split_tibble(tfbs, "sequence_name")

  for (i in 1:nrow(human_retro)) {
    retro_id <- as.character(human_retro[i, "retroID"])
    pg_id <- as.character(human_retro[i, "parentID"])

    retro_motif <- s[retro_id][[1]][, "motif_id"]
    pg_motif <- s[pg_id][[1]][, "motif_id"]

    if (length(intersect(retro_motif$motif_id, pg_motif$motif_id)) == 0) {
      human_retro[i, "n_shared_motif"] <- 0
      human_retro[i, "shared_motif"] <- ""
    } else {
      shared <- vintersect(retro_motif$motif_id, pg_motif$motif_id) # vintersect with duplicates
      human_retro[i, "n_shared_motif"] <- length(shared)
      uniq_motif <-
        human_retro[i, "shared_motif"] <- paste(unlist(as.list(shared)), collapse = ", ")
    }
  }

  df <- human_retro %>% select(c(retroID, GM12878, HMEC, HUVEC, IMR90, K562, KBM7, NHEK))
  df <- melt(df, id.vars = "retroID", variable.name = "cell_line", value.name = "colocal")
  df <- df %>% left_join(human_retro, by = "retroID")

  ggplot(df, aes(x = colocal, y = n_shared_motif)) +
    geom_violin(aes(fill = factor(colocal))) +
    geom_boxplot(width = 0.3, outlier.shape = NA, size = 0.8) +
    facet_grid(. ~ cell_line) +
    geom_signif(
      comparisons = list(c("Colocalized", "Noncolocalized")),
      map_signif_level = TRUE, size = 0.6, textsize = 8, vjust = 0.5
    ) +
    labs(x = "", y = "Number of shared TFBSs") +
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
shared_tfbs_plot("data/TFBSs/fimo_TSS_3kb_q0.01.tsv")
ggsave("plots/fig3C.pdf", width = 12, height = 8)

###
p1 <- shared_tfbs_plot("data/TFBSs/fimo_TSS_1kb_q0.01.tsv")
p2 <- shared_tfbs_plot("data/TFBSs/fimo_TSS_up3kb_q0.01.tsv")

library(ggpubr)
ggarrange(p1, p2, common.legend = TRUE)
ggsave("plots/figS15.pdf", width = 20, height = 10)



##### Compare retrocopies and parental genes #####

tfbs_plot <- function(tfbs, type = c("Retrocopies", "Parental genes")) {
  if (type == "Retrocopies") {
    retro <- tfbs %>%
      filter(grepl("hsa_", sequence_name)) %>%
      arrange(sequence_name)
    retro_tfbs_count <- as.data.frame(table(retro$sequence_name))
    colnames(retro_tfbs_count) <- c("retroID", "TFBS_count")

    df <- human_retro %>% select(c(retroID, GM12878, HMEC, HUVEC, IMR90, K562, KBM7, NHEK))
    df <- melt(df, id.vars = "retroID", variable.name = "cell_line", value.name = "colocal") %>%
      left_join(retro_tfbs_count, by = "retroID")
  } else if (type == "Parental genes") {
    pg <- tfbs %>%
      filter(grepl("ENSG", sequence_name)) %>%
      arrange(sequence_name)
    pg_tfbs_count <- as.data.frame(table(pg$sequence_name))
    colnames(pg_tfbs_count) <- c("parentID", "TFBS_count")

    df <- human_retro %>% select(c(parentID, GM12878, HMEC, HUVEC, IMR90, K562, KBM7, NHEK))
    df <- melt(df, id.vars = "parentID", variable.name = "cell_line", value.name = "colocal") %>%
      left_join(pg_tfbs_count, by = "parentID")
  }
  ggplot(df, aes(x = colocal, y = log2(TFBS_count))) +
    geom_violin(aes(fill = factor(colocal))) +
    geom_boxplot(width = 0.3, outlier.shape = NA, size = 0.8) +
    facet_grid(. ~ cell_line) +
    geom_signif(
      comparisons = list(c("Colocalized", "Noncolocalized")),
      map_signif_level = TRUE, size = 0.6, textsize = 8, vjust = 0.5
    ) +
    labs(x = "", y = bquote(bold(log[2] * "(#TFBSs)")), tag = "", title = type) +
    guides(fill = guide_legend(title = NULL)) +
    theme_publication(baseSize = 24) +
    theme(
      legend.position = "top",
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    ) +
    scale_fill_discrete(labels = c("Colocalized", "Noncolocalized"))
}

tfbs <- read_tsv("data/TFBSs/fimo_TSS_3kb_q0.01.tsv")

p1 <- tfbs_plot(tfbs, "Retrocopies")
p2 <- tfbs_plot(tfbs, "Parental genes")

library(ggpubr)
ggarrange(p1, p2,
  labels = "AUTO", common.legend = TRUE, legend = "bottom",
  font.label = list(size = 28, color = "black", face = "bold", family = NULL)
)
ggsave("plots/figS14.pdf", width = 20, height = 10)



##### GTEx data analysis #####

### Read GTEx data
library(tidyverse)
load("Rdata/2-colocalization_250kb.Rdata")

human_retro$colocNum <- apply(human_retro, 1, function(x) length(which(x == "Colocalized")))
human_retro <- human_retro %>%
  dplyr::mutate(colocType = ifelse(colocNum >= 2, "≥2", as.character(colocNum)))

# gene
gtex_gene <- read_tsv("data/GTEx/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct",
  skip = 2
)
gtex_gene <- gtex_gene %>%
  dplyr::rename(parentID = Name) %>%
  dplyr::mutate(parentID = substr(parentID, 1, 15))

human_retro <- human_retro %>% left_join(gtex_gene, by = "parentID")


# trx (needs a lot of memory)
gtex_trx <- read_tsv("data/GTEx/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct",
  skip = 2
)
gtex_trx <- gtex_trx %>%
  dplyr::rename(annotatedRetroID = transcript_id) %>%
  dplyr::mutate(annotatedRetroID = substr(annotatedRetroID, 1, 15))

human_retro <- human_retro %>% left_join(gtex_trx, by = "annotatedRetroID")

rm(gtex_trx, gtex_gene)
save.image("Rdata/6-GTEx.Rdata")


###
human_retro <- human_retro %>%
  dplyr::filter(!is.na(Description), !is.na(gene_id))

library(doParallel)
# Setup backend to use many processors
totalCores <- detectCores()

# Leave one core to avoid overload your computer
cluster <- makeCluster(totalCores[1] - 1)
registerDoParallel(cluster)

coor <- foreach(i = 1:3635, .combine = rbind) %dopar% {
  x <- as.numeric(human_retro[i, 24:17405])
  y <- as.numeric(human_retro[i, 17407:34788])
  corr <- c(corr, cor.test(x, y)$estimate)
}

corr_df <- as.data.frame(corr)
df <- human_retro %>% dplyr::select(retroID:colocType)
corr_df <- cbind(df, corr_df)

corr_df$colocType <- factor(corr_df$colocType, levels = c("0", "1", "≥2"))

source("theme_publication.R")
ggplot(corr_df, aes(x = corr)) +
  stat_ecdf(aes(color = as.factor(colocType)),
    geom = "step", linewidth = 1.5
  ) +
  labs(
    x = bquote(bold("Correlation coefficient")), y = "Cumulative fraction"
  ) +
  annotate("text",
    x = -0.1, y = 0.9,
    label = "0 vs 1: p = 2.2e-02\n1 vs 2: p = 3.5e-01\n0 vs 2: p = 8.3e-04", size = 6
  ) +
  theme_publication(baseSize = 24) +
  theme(legend.position = "top") +
  guides(color = guide_legend(title = "No. of colocalized cell lines"))

ggsave("plots/figS16.pdf", width = 8, height = 8, device = cairo_pdf)
