##### ENCODE RNA Dashboard data analysis #####

##### Clean data #####
load("Rdata/2-colocalization_250kb.Rdata")
library(tidyverse)

cell_lines <- c("GM12878", "HUVEC", "K562", "NHEK")
folder <- "data/ENCODE_RNA_Dashboard/"

rna_results <- vector()
for (cell_line in cell_lines) {
  rna_results <- append(rna_results, list.files(
    path = paste0(folder, cell_line),
    recursive = FALSE, full.names = TRUE
  ))
}


# read each sample and append to the human_retro dataframe
for (f in rna_results) {
  base_name <- str_match(f, "RnaSeq\\s*(.*?)\\s*_v19.tsv")[, 2]
  rna <- read_tsv(f)

  rna_R <- rna %>%
    select("Gene Name", FPKM) %>%
    rename(retroID = "Gene Name")

  rna_P <- rna %>%
    select("Gene ID", FPKM) %>%
    rename(parentID = "Gene ID")
  rna_P$parentID <- str_sub(rna_P$parentID, 1, 15)

  human_retro <- human_retro %>%
    left_join(rna_R, by = "retroID") %>%
    left_join(rna_P, by = "parentID") %>%
    rename(!!paste(base_name, "_FPKM.R", sep = "") := "FPKM.x") %>%
    rename(!!paste(base_name, "_FPKM.P", sep = "") := "FPKM.y")
}


# get the mean FPKM of different replicates
for (cell_line in cell_lines) {
  cell_line <- str_to_sentence(cell_line)
  for (compartment in c("Cell", "Nucleus", "Cytosol")) {
    for (tech in c("Pap", "Longnonpolya")) {
      for (gene in c("R", "P")) {
        human_retro <- human_retro %>%
          mutate(!!paste(cell_line, compartment, tech, ".", gene, sep = "") :=
            select(
              human_retro,
              matches(paste(cell_line, compartment, tech, ".*FPKM.", gene, sep = ""))
            ) %>%
            rowMeans())
      }
    }
  }
}


human_retro <- human_retro %>%
  select(!matches("AlnRep")) %>%
  select_if(function(x) {
    !all(is.na(x))
  })

save(human_retro, cell_lines, file = "Rdata/5-expression.RData")



##### Plots #####

load("Rdata/5-expression.RData")
source("theme_publication.R")
library(reshape2)
library(tidyverse)
library(ggsignif)

###
ecdf_plot <- function(cell_line, component, title, text_x, text_y) {
  x <- paste0(cell_line, component)

  C <- human_retro %>% filter(!!as.symbol(toupper(cell_line)) == "Colocalized")
  N <- human_retro %>% filter(!!as.symbol(toupper(cell_line)) == "Noncolocalized")

  t <- wilcox.test(log10(C[, x][[1]] + 1), log10(N[, x][[1]] + 1))

  ggplot(human_retro, aes_string(x = paste0("log10(", x, "+ 1)"))) +
    stat_ecdf(aes_string(color = toupper(cell_line)),
      geom = "step", size = 1.5
    ) +
    scale_color_manual(
      values = c(rgb(231 / 255, 126 / 255, 114 / 255), rgb(85 / 255, 188 / 255, 194 / 255)),
      name = "", labels = c("Colocalized", "Noncolocalized")
    ) +
    annotate("text", x = text_x, y = text_y, label = sprintf("p = %.1e", t$p.value), size = 8) +
    labs(
      x = bquote(bold(log[10] * "(FPKM + 1)")), y = "Cumulative fraction",
      tag = "", title = title
    ) +
    theme_publication(baseSize = 24)
}


p1 <- ecdf_plot("Gm12878", "CellPap.R", "Retrocopies (cell ployA+)", 1.5, 0.7)
p2 <- ecdf_plot("Gm12878", "CellPap.P", "Parental genes (cell ployA+)", 1, 1)
p3 <- ecdf_plot("Gm12878", "CellLongnonpolya.R", "Retrocopies (cell ployA-)", 1.5, 0.7)
p4 <- ecdf_plot("Gm12878", "CellLongnonpolya.P", "Parental genes (cell ployA-)", 1, 1)


library(ggpubr)
ggarrange(p1, p2, p3, p4, common.legend = TRUE)
ggsave("plots/fig3A.pdf", width = 11, height = 11)


###
exp_plot <- function(component, gene, y_pos, y_lim, tip, title) {
  # clean data
  for (cell_line in cell_lines) {
    cell_line <- str_to_sentence(cell_line)
    if (gene == "retro") {
      fpkm <- paste0(cell_line, component, ".R", sep = "")
      df <- human_retro %>%
        select(c("retroID", toupper(cell_line))) %>%
        melt(id.vars = "retroID", variable.name = "cell_line", value.name = "coloc") %>%
        left_join(human_retro[, c("retroID", fpkm)], by = "retroID") %>%
        rename(FPKM = fpkm)
    } else if (gene == "parent") {
      fpkm <- paste0(cell_line, component, ".P", sep = "")
      df <- human_retro %>%
        select(c("parentID", toupper(cell_line))) %>%
        melt(id.vars = "parentID", variable.name = "cell_line", value.name = "coloc") %>%
        left_join(human_retro[, c("parentID", fpkm)], by = "parentID") %>%
        rename(FPKM = fpkm)
      df <- distinct(df)
    }
    assign(cell_line, df)
  }

  df <- rbind(Gm12878, Huvec, K562, Nhek)

  # plot
  ggplot(df, aes(x = coloc, y = log2(FPKM + 1))) +
    geom_boxplot(aes(fill = factor(coloc)), outlier.shape = NA, size = 1.1) +
    facet_grid(. ~ cell_line) +
    geom_signif(
      comparisons = list(c("Colocalized", "Noncolocalized")), y_position = c(y_pos),
      tip_length = c(tip, tip),
      map_signif_level = TRUE, size = 0.8, textsize = 10
    ) +
    labs(x = "", y = bquote(bold(log[2] * "(FPKM + 1)")), tag = "", title = title) +
    ylim(0, y_lim) +
    guides(fill = guide_legend(title = NULL)) +
    theme_publication(baseSize = 32) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    ) +
    scale_fill_discrete(labels = c("Colocalized", "Noncolocalized"))
}


# Retrocopies
p1 <- exp_plot("CellPap", "retro", 2.2, 3, 0.008, "Cell polyA+")
p2 <- exp_plot("NucleusPap", "retro", 1.2, 1.8, 0.008, "Nucleus polyA+")
p3 <- exp_plot("CytosolPap", "retro", 1.1, 1.8, 0.006, "Cytosol polyA+")
p4 <- exp_plot("CellLongnonpolya", "retro", 2.2, 3, 0.01, "Cell polyA-")
p5 <- exp_plot("NucleusLongnonpolya", "retro", 1.5, 2, 0.01, "Nucleus polyA-")
p6 <- exp_plot("CytosolLongnonpolya", "retro", 0.35, 0.8, 0.004, "Cytosol polyA-")

library(ggpubr)
ggarrange(p1, p2, p3, p4, p5, p6, nrow = 2, ncol = 3, common.legend = TRUE)
ggsave("plots/figS12.pdf", width = 28, height = 20)

# Parental genes
p1 <- exp_plot("CellPap", "parent", 13, 15, 0.035, "Cell polyA+")
p2 <- exp_plot("NucleusPap", "parent", 13, 15, 0.05, "Nucleus polyA+")
p3 <- exp_plot("CytosolPap", "parent", 13, 15, 0.05, "Cytosol polyA+")
p4 <- exp_plot("CellLongnonpolya", "parent", 13, 15, 0.04, "Cell polyA-")
p5 <- exp_plot("NucleusLongnonpolya", "parent", 9, 10, 0.04, "Nucleus polyA-")
p6 <- exp_plot("CytosolLongnonpolya", "parent", 8.5, 10, 0.04, "Cytosol polyA-")

ggarrange(p1, p2, p3, p4, p5, p6, nrow = 2, ncol = 3, common.legend = TRUE)
ggsave("plots/figS13.pdf", width = 28, height = 20)


###
exp_plot2 <- function(fraction, title) {
  for (cell_line in cell_lines) {
    cell_line <- str_to_sentence(cell_line)
    rna <- paste0(cell_line, fraction, ".R")

    df <- human_retro %>%
      mutate(exp_class := ifelse(!!sym(rna) >= 1, "FPKM ≥ 1",
        ifelse(!!sym(rna) < 0.1, "FPKM < 0.1", "0.1 ≤ FPKM < 1")
      )) %>%
      count(exp_class) %>%
      mutate(pct = prop.table(n) * 100) %>%
      drop_na() %>%
      cbind(data.frame(cell_line = rep(toupper(cell_line), 3)))

    assign(cell_line, df)
  }

  df <- rbind(Gm12878, Huvec, K562, Nhek)
  df$exp_class <- factor(df$exp_class, levels = c("FPKM ≥ 1", "0.1 ≤ FPKM < 1", "FPKM < 0.1"))

  # plot
  ggplot(df, aes(x = cell_line, y = pct, fill = exp_class)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = paste0(sprintf("%1.1f", pct), "%", "\n(", n, ")")),
      position = position_stack(vjust = 0.5), size = 6
    ) +
    labs(x = "", y = "Percentage", title = title) +
    theme_publication(baseSize = 28) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    theme(axis.title = element_text(size = 22)) +
    scale_fill_discrete(name = "") +
    theme(legend.position = "top")
}

p1 <- exp_plot2("CellPap", "Cell polyA+")
p2 <- exp_plot2("NucleusPap", "Nucleus polyA+")
p3 <- exp_plot2("CytosolPap", "Cytosol polyA+")
p4 <- exp_plot2("CellLongnonpolya", "Cell polyA-")
p5 <- exp_plot2("NucleusLongnonpolya", "Nucleus polyA-")
p6 <- exp_plot2("CytosolLongnonpolya", "Cytosol polyA-")

library(ggpubr)
ggarrange(p1, p2, p3, p4, p5, p6, nrow = 2, ncol = 3, common.legend = TRUE)
ggsave("plots/figS11.pdf", width = 20, height = 16, device = cairo_pdf)
