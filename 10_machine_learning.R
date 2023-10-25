##### Logistic regression #####

library(tidyverse)

for (cell_line in c("GM12878", "HUVEC", "K562", "NHEK")) {
  assign(
    paste0(cell_line, "_df"),
    read_tsv(paste0("data/Machine_Learning/", cell_line, "_features.tab"))
  )
}


# logistic regression
library(caTools)
library(caret)

# GM12878_df$GM12878 <- as.factor(GM12878_df$GM12878)

set.seed(100)
sample_split <- sample.split(Y = GM12878_df$GM12878, SplitRatio = 0.7)
train_set <- subset(x = GM12878_df, sample_split == TRUE)
test_set <- subset(x = GM12878_df, sample_split == FALSE)

logistic <- glm(GM12878 ~ ., data = train_set, family = "binomial")

library(pROC)
test_prob <- predict(logistic, newdata = test_set, type = "response")
par(pty = "s")
test_roc <- roc(test_set$GM12878 ~ test_prob, plot = TRUE, print.auc = TRUE, legacy.axes = TRUE)
# Plot size: 6 x 6; AUC: 0.785


# Importance
library(Boruta)
boruta_output <- Boruta(GM12878 ~ ., data = train_set, doTrace = 0)

par(mar = c(12, 2, 2, 2))
plot(boruta_output, ces.axis = 0.5, las = 2, xlab = "", main = "Feature importance")
# Plot size: 12 x 12
