# download and modified from https://rpubs.com/Koundy/71792
# by Koundinya Desiraju

theme_publication <- function(baseSize = 16) {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size = baseSize)
  + theme(
      plot.title = element_text(
        face = "bold",
        size = rel(1), hjust = 0.5
      ),
      text = element_text(),
      panel.background = element_rect(colour = NA),
      plot.background = element_rect(colour = NA),
      panel.border = element_rect(colour = NA),
      axis.title = element_text(face = "bold", size = rel(1)),
      axis.title.y = element_text(angle = 90, vjust = 2),
      axis.title.x = element_text(vjust = -0.2),
      axis.text = element_text(),
      axis.line = element_line(),
      axis.ticks = element_line(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.key = element_rect(colour = NA),
      legend.key.size = unit(0.5, "cm"),
      plot.margin = unit(c(10, 5, 5, 5), "mm"),
      strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
      strip.text = element_text(face = "bold")
    ))
}
