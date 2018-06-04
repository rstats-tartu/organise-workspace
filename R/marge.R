
#' check and install required packages
install.packages("pacman")
pacman::p_load(tidyverse, readxl)

#' load required libraries
library(tidyverse)
library(readxl)

#' import data file from excel file
#' get sheet names
(sheets <- excel_sheets("docs/Statistika_Küsimused.xlsx"))
#' import region of interest from sheet
marge <- read_excel("docs/Statistika_Küsimused.xlsx", 
                     sheet = "Marge",
                     range = c("C16:L31"))

#' we need to fix colnames
colnames(marge) <- c("chem", "ec50", "algae", "water", "ph", "org", "cond", "nano_size", "nano_z", "solubility")

#' fix typos
marge <- mutate_if(marge, is.character, str_trim) %>% 
  mutate(water = case_when(
    str_detect(water, "Ülemsite") ~ "Ülemiste",
               TRUE ~ water
  ))

#' we can test the importance of each predictor variable by using linear model
#' But first we need to throw out columns with NA-s
#' another problem is that every row is unique, meaning that we don't have
#' enough replicates. We have to leave it as it is and collect more data
marge_sub <- dplyr::select(marge, -starts_with("nano"), -cond)
m1 <- lm(ec50 ~ ., data = marge_sub)
summary(m1)

library(sandwich)
library(lmtest)
coeftest(m1, vcov = sandwich)

library(multcomp)
summary(glht(m1))

#' PCA can be used as a means of extracting dominant patterns in a group of predictor
#' variables in the later building of a model, basically you reduce the dimensions of
#' your dataset
marge_num <- mutate_if(marge_sub, is.character, function(x) as.numeric(as.factor(x)))
pca <- prcomp(marge_num, retx = TRUE, center = TRUE, scale = TRUE)

expl.var <- round(pca$sdev^2 / sum(pca$sdev^2) * 100) # percent explained variance

pca#' plot first two principal components
plot(pca$x[,c(1, 2)], col = c(1:4)[marge_num$chem], pch = 16,
     xlab = paste0("PC2 (", expl.var[1], "%)"), 
     ylab = paste0("PC1 (", expl.var[2], "%)")) 
legend("topright", 
       legend = unique(marge_sub$chem), 
       fill = c(2:4), 
       border = c(2:4))
