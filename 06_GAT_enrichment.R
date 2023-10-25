##### GAT enrichment analysis #####

library(tidyverse)
library(reshape2)
library(ggpubr)
library(formattable)
source("theme_publication.R")


##### Plot fold changes for Hi-C annotation in 5 cell lines #####

# Features analyzed using GAT
features <- c(
  "Annotation", "Protein-coding genes",
  "Unprocessed pseudogenes", "Parental genes", "Retrocopies",
  "Colocalized parental genes", "Noncolocalized parental genes", "retroCNVs",
  "retroCNVs colocalized in 0 cell line", "retroCNVs colocalized in 1 cell line",
  "retroCNVs colocalized in ≥ 2 cell lines",
  "Colocalized retrocopies", "Intergenic colocalized retrocopies",
  "Intragenic colocalized retrocopies",
  "Noncolocalized retrocopies", "Intergenic noncolocalized retrocopies",
  "Intragenic noncolocalized retrocopies"
)

factors <- c(
  "Retrocopies", "Colocalized retrocopies", "Intergenic colocalized retrocopies",
  "Intragenic colocalized retrocopies", "Noncolocalized retrocopies",
  "Intergenic noncolocalized retrocopies", "Intragenic noncolocalized retrocopies",
  "Parental genes", "Colocalized parental genes",
  "Noncolocalized parental genes", "Protein-coding genes", "Unprocessed pseudogenes",
  "retroCNVs", "retroCNVs colocalized in 0 cell line", "retroCNVs colocalized in 1 cell line",
  "retroCNVs colocalized in ≥ 2 cell lines"
)


for (cell_line in c("GM12878", "HMEC", "HUVEC", "IMR90", "K562", "K562_SPIN")) {
  # fold changes were computed using GAT
  fc <- read_tsv(paste0("data/GAT_enrichment/", cell_line, "_fold_enrichment.txt"))
  colnames(fc) <- features
  fc <- melt(fc, id.vars = "Annotation", variable.name = "Type", value.name = "Fold change")
  fc$Type <- factor(fc$Type, levels = rev(factors))
  fc$`Fold change` <- round(fc$`Fold change`, digits = 2)

  if (cell_line == "K562_SPIN") {
    fc_max <- 3.3
    angle <- 45
    hjust <- 1
  } else {
    fc_max <- max(fc$`Fold change`)
    angle <- 0
    hjust <- 0.5
  }

  p <- ggplot(fc, aes(Annotation, Type)) +
    geom_tile(aes(fill = `Fold change`), color = "black") +
    geom_text(aes(label = formattable::digits(`Fold change`, digits = 2)), parse = FALSE) +
    scale_fill_gradient2(
      low = "blue", mid = "white", high = "red",
      midpoint = 1, limits = c(0, fc_max)
    ) +
    labs(x = "", y = "", title = cell_line) +
    theme_publication() +
    theme(
      axis.ticks.x = element_line(colour = NA),
      axis.ticks.y = element_line(colour = NA),
      axis.line = element_line(colour = NA),
      axis.text.x = element_text(angle = angle, vjust = 1, hjust = hjust)
    )

  assign(paste0("p_", cell_line), p)
}


# Save plot for GM12878
ggsave("plots/fig2A.pdf", plot = p_GM12878, width = 9.7, height = 6, device = cairo_pdf)

# Save plot for HMEC, HUVEC, IMR90, and K562
ggarrange(p_HMEC, p_HUVEC, p_IMR90, p_K562,
  labels = "AUTO",
  font.label = list(size = 24, face = "bold")
)
ggsave("plots/figS6.pdf", width = 18, height = 10, device = cairo_pdf)

# Save plot for K562 SPIN enrichment
ggsave("plots/figS9.pdf", plot = p_K562_SPIN, width = 9.7, height = 6, device = cairo_pdf)



##### Sankey plots #####
load("Rdata/2-colocalization_250kb.Rdata")
source("theme_publication.R")

sankey_plot <- function(df, title) {
  require(ggsankey)

  factors <- c("A1", "A2", "B1", "B2", "B3", "NA")
  df$node <- factor(df$node, levels = rev(factors))
  df$next_node <- factor(df$next_node, levels = rev(factors))

  ggplot(df, aes(
    x = x,
    next_x = next_x,
    node = node,
    next_node = next_node,
    fill = factor(node),
    label = node
  )) +
    geom_sankey(
      flow.alpha = 0.5,
      node.color = "black",
      show.legend = TRUE
    ) +
    geom_sankey_label(size = 4, color = "black", fill = "white", hjust = -0.4) +
    scale_fill_manual(values = c(
      "NA" = "gray",
      "A1" = "deeppink1", "A2" = "salmon",
      "B1" = "springgreen3", "B2" = "darkturquoise", "B3" = "dodgerblue"
    )) +
    labs(title = title) +
    theme_publication() +
    theme(
      axis.title.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      axis.line = element_blank(),
      legend.title = element_blank(),
      legend.position = "none"
    )
}


sankey_plotRP <- function(R, P, Cell) {
  require(tidyverse)
  require(ggsankey)
  require(ggpubr)

  compR <- read_tsv(R, col_names = F)
  names(compR) <- c("retroID", "Retrocopy")

  compP <- read_tsv(P, col_names = F)
  names(compP) <- c("parentID", "Parental gene")

  # Colocalized pairs
  dfC <- human_retro %>%
    filter(!!sym(Cell) == "Colocalized") %>%
    left_join(compR, by = "retroID") %>%
    left_join(compP, by = "parentID")

  dfC <- dfC %>%
    select(Retrocopy, `Parental gene`) %>%
    replace_na(list(Retrocopy = "NA", `Parental gene` = "NA"))

  dfC <- dfC %>% make_long(Retrocopy, `Parental gene`)

  p1 <- sankey_plot(dfC, "Colocalized")


  # Noncolocalized pairs
  dfN <- human_retro %>%
    filter(!!sym(Cell) == "Noncolocalized") %>%
    left_join(compR, by = "retroID") %>%
    left_join(compP, by = "parentID")

  dfN <- dfN %>%
    select(Retrocopy, `Parental gene`) %>%
    replace_na(list(Retrocopy = "NA", `Parental gene` = "NA"))

  dfN <- dfN %>% make_long(Retrocopy, `Parental gene`)

  p2 <- sankey_plot(dfN, "Noncolocalized")

  # save
  return(ggarrange(p1, p2))
}


for (cell_line in c("GM12878", "HMEC", "HUVEC", "IMR90", "K562")) {
  p <- sankey_plotRP(
    paste0("data/GAT_enrichment/", cell_line, "_subcompartments_retro.txt"),
    paste0("data/GAT_enrichment/", cell_line, "_subcompartments_parent.txt"),
    cell_line
  )

  assign(paste0("sankey_", cell_line), p)
}


# Save plot for GM12878
ggsave("plots/fig2B.pdf", plot = sankey_GM12878, width = 9.7, height = 6, device = cairo_pdf)

# Save plot for HMEC, HUVEC, IMR90, and K562
ggarrange(sankey_HMEC, sankey_HUVEC, sankey_IMR90, sankey_K562,
  labels = "AUTO",
  font.label = list(size = 24, face = "bold")
)
ggsave("plots/figS7.pdf", width = 20, height = 12, device = cairo_pdf)
