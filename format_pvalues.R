library(ggplot2)
library(ggsignif)

# A function to format a numeric p-value in scientific notation
# max_pd parameter refers to max decimal places before switching to scientific notation.
# Otherwise the full decimal is displayed.
# p-values are given to one significant digit

pform <- function(p, max_dp = 3) {
  # map the superscript versions of numbers
  # (obtained from https://en.wikipedia.org/wiki/Unicode_subscripts_and_superscripts)
  ucmap <- c(
    "0" = "\u2070",
    "1" = "\u00B9",
    "2" = "\u00B2",
    "3" = "\u00B3",
    "4" = "\u2074",
    "5" = "\u2075",
    "6" = "\u2076",
    "7" = "\u2077",
    "8" = "\u2078",
    "9" = "\u2079",
    "-" = "\u207B"
  )

  pe <- formatC(signif(p, 2), digits = 1, format = "e")

  pv <- unlist(strsplit(pe, split = "e"))

  if (abs(as.numeric(pv[2])) > max_dp) {
    pv2 <- as.character(as.numeric(pv[2]))

    pv2v <- unlist(strsplit(pv2, split = ""))

    ps <- paste(ucmap[pv2v], sep = "", collapse = "")

    paste(pv[1], "\u00D7", "10", ps, collapse = "", sep = "")
  } else {
    formatC(signif(p, 1), digits = max_dp, format = "f", drop0trailing = T)
  }
}
