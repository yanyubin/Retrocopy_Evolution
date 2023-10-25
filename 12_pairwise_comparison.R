##### Analysis twin retrocopies with one colocalized and the other noncolocalized #####

library(tidyverse)
library(reshape2)
library(ggpubr)
library(rstatix)

twin_retros <- "data/Conservation_Score/parent_with_2retro.tsv"

retro <- read_tsv(twin_retros)


##### expression
C <- read_tsv(twin_retros) %>%
  filter(GM12878 == "C") %>%
  select(retroID, parentID, GM12878, Gm12878CellTotal.R) %>%
  mutate(GM12878 = "Colocalized")


N <- read_tsv(twin_retros) %>%
  filter(GM12878 == "N") %>%
  select(retroID, parentID, GM12878, Gm12878CellTotal.R) %>%
  mutate(GM12878 = "Noncolocalized")

d <- N %>%
  left_join(C, by = "parentID") %>%
  filter(!is.na(retroID.y))

d1 <- d %>% select(retroID.x, GM12878.x, Gm12878CellTotal.R.x)
names(d1) <- c("retroID", "GM12878", "value")

d2 <- d %>% select(retroID.y, GM12878.y, Gm12878CellTotal.R.y)
names(d2) <- c("retroID", "GM12878", "value")


df <- rbind(d1, d2) %>% mutate(value = log2(value + 1))

df$GM12878 <- factor(df$GM12878, levels = c("Colocalized", "Noncolocalized"), ordered = TRUE)


stat.test <- df %>%
  t_test(value ~ GM12878, paired = TRUE) %>%
  add_significance()

stat.test <- stat.test %>% add_xy_position(x = "GM12878")


p1 <- ggpaired(df,
  x = "GM12878", y = "value",
  fill = c(rgb(231 / 255, 126 / 255, 114 / 255), rgb(85 / 255, 188 / 255, 194 / 255)),
  line.color = "gray", line.size = 0.4
) +
  labs(x = NULL, y = bquote(bold(log[2] * "(FPKM + 1)")), title = "RNA-seq") +
  stat_pvalue_manual(stat.test, label = "{p.signif}", size = 6) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.10))) +
  theme(text = element_text(size = 16)) +
  theme(plot.title = element_text(hjust = 0.5))



##### no. of shared motifs
C <- read_tsv(twin_retros) %>%
  filter(GM12878 == "C") %>%
  select(retroID, parentID, GM12878, numSharedMotifs) %>%
  mutate(GM12878 = "Colocalized")


N <- read_tsv(twin_retros) %>%
  filter(GM12878 == "N") %>%
  select(retroID, parentID, GM12878, numSharedMotifs) %>%
  mutate(GM12878 = "Noncolocalized")

d <- N %>%
  left_join(C, by = "parentID") %>%
  filter(!is.na(retroID.y))

d1 <- d %>% select(retroID.x, GM12878.x, numSharedMotifs.x)
names(d1) <- c("retroID", "GM12878", "value")

d2 <- d %>% select(retroID.y, GM12878.y, numSharedMotifs.y)
names(d2) <- c("retroID", "GM12878", "value")


df <- rbind(d1, d2)

df$GM12878 <- factor(df$GM12878, levels = c("Colocalized", "Noncolocalized"), ordered = TRUE)


stat.test <- df %>%
  t_test(value ~ GM12878, paired = TRUE) %>%
  add_significance()

stat.test <- stat.test %>% add_xy_position(x = "GM12878")


p2 <- ggpaired(df,
  x = "GM12878", y = "value",
  fill = c(rgb(231 / 255, 126 / 255, 114 / 255), rgb(85 / 255, 188 / 255, 194 / 255)),
  line.color = "gray", line.size = 0.4
) +
  labs(x = NULL, y = "Number of shared TFBSs", title = "Shared regulatory elements") +
  stat_pvalue_manual(stat.test, label = "{p.signif}", size = 6) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.10))) +
  theme(text = element_text(size = 16)) +
  theme(plot.title = element_text(hjust = 0.5))



##### Ribo-seq
C <- read_tsv(twin_retros) %>%
  filter(GM12878 == "C") %>%
  select(retroID, parentID, GM12878, maxRiboRPKM) %>%
  mutate(GM12878 = "Colocalized") %>%
  filter(maxRiboRPKM >= 1)


N <- read_tsv(twin_retros) %>%
  filter(GM12878 == "N") %>%
  select(retroID, parentID, GM12878, maxRiboRPKM) %>%
  mutate(GM12878 = "Noncolocalized") %>%
  filter(maxRiboRPKM >= 1)

d <- N %>%
  left_join(C, by = "parentID") %>%
  filter(!is.na(retroID.y))

d1 <- d %>% select(retroID.x, GM12878.x, maxRiboRPKM.x)
names(d1) <- c("retroID", "GM12878", "value")

d2 <- d %>% select(retroID.y, GM12878.y, maxRiboRPKM.y)
names(d2) <- c("retroID", "GM12878", "value")


df <- rbind(d1, d2) %>% mutate(value = log10(value))

df$GM12878 <- factor(df$GM12878, levels = c("Colocalized", "Noncolocalized"), ordered = TRUE)


stat.test <- df %>%
  t_test(value ~ GM12878, paired = TRUE) %>%
  add_significance()

stat.test <- stat.test %>% add_xy_position(x = "GM12878")


p3 <- ggpaired(df,
  x = "GM12878", y = "value",
  fill = c(rgb(231 / 255, 126 / 255, 114 / 255), rgb(85 / 255, 188 / 255, 194 / 255)),
  line.color = "gray", line.size = 0.4
) +
  labs(x = NULL, y = bquote(bold(log[10] * "(RPKM)")), title = "Ribo-seq") +
  stat_pvalue_manual(stat.test, label = "{p.signif}", size = 6) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.10))) +
  theme(text = element_text(size = 16)) +
  theme(plot.title = element_text(hjust = 0.5))



##### phastCons
C <- read_tsv(twin_retros) %>%
  filter(GM12878 == "C") %>%
  select(retroID, parentID, GM12878, phastConsScore) %>%
  mutate(GM12878 = "Colocalized")


N <- read_tsv(twin_retros) %>%
  filter(GM12878 == "N") %>%
  select(retroID, parentID, GM12878, phastConsScore) %>%
  mutate(GM12878 = "Noncolocalized")

d <- N %>%
  left_join(C, by = "parentID") %>%
  filter(!is.na(retroID.y))

d1 <- d %>% select(retroID.x, GM12878.x, phastConsScore.x)
names(d1) <- c("retroID", "GM12878", "value")

d2 <- d %>% select(retroID.y, GM12878.y, phastConsScore.y)
names(d2) <- c("retroID", "GM12878", "value")


df <- rbind(d1, d2)
df$GM12878 <- factor(df$GM12878, levels = c("Colocalized", "Noncolocalized"), ordered = TRUE)


stat.test <- df %>%
  t_test(value ~ GM12878, paired = TRUE) %>%
  add_significance()

stat.test <- stat.test %>% add_xy_position(x = "GM12878")


p4 <- ggpaired(df,
  x = "GM12878", y = "value",
  fill = c(rgb(231 / 255, 126 / 255, 114 / 255), rgb(85 / 255, 188 / 255, 194 / 255)),
  line.color = "gray", line.size = 0.4
) +
  labs(x = NULL, y = "PhastCons score", title = "Conservation") +
  stat_pvalue_manual(stat.test, label = "{p.signif}", size = 6) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.10))) +
  theme(text = element_text(size = 16)) +
  theme(plot.title = element_text(hjust = 0.5))



library(ggpubr)
ggarrange(p1, p2, p3, p4,
  labels = "AUTO",
  font.label = list(size = 18, color = "black", face = "bold", family = "sans")
)

ggsave("plots/figS30.pdf", width = 10, height = 10)
