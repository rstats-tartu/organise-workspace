
#' check and install required packages
# install.packages("pacman")
pacman::p_load(tidyverse, readxl, car, AER, MASS, modelr)

#' load required libraries
library(tidyverse)
library(readxl)

#' import data file from excel file
#' get sheet names
(sheets <- excel_sheets("docs/Statistika_Küsimused.xlsx"))
#' import region of interest from sheet
margit <- read_excel("docs/Statistika_Küsimused.xlsx", 
                       sheet = "Margit",
                       range = c("B5:Z50"))

# dates <- margit[1,]
# dates <- unlist(dates)
#' we need conc and number of pedigree and length
msml <- margit[5:nrow(margit), c(1, 2, 24, 25)]

#' remove rows with only NAs and convert data to numeric
msml_filt <- filter(msml, apply(msml, 1, function(x) !all(is.na(x))))
colnames(msml_filt) <- c("conc", "id", "n", "len")  
msml_filt <- mutate(msml_filt, 
                    conc = case_when(
                      str_detect(conc, "\\(konts") ~ "NA",
                      conc == "KONTROLL" ~ "1",
                      TRUE ~ conc)) %>%
  mutate_at(vars("conc", "n", "len"), as.numeric) %>% 
  fill(conc)

#' let's plot this dataset
ggplot(data = msml_filt) +
  geom_point(mapping = aes(x = log10(conc), y = n)) +
  stat_summary(fun.data = mean_se, mapping = aes(x = log10(conc), y = n), color = "red") +
  labs(x = "Particles/ml",
       y = "Number of pedigree")

ggplot(data = msml_filt) +
  geom_point(mapping = aes(x = log10(conc), y = len)) +
  stat_summary(fun.data = mean_se, mapping = aes(x = log10(conc), y = len), color = "red") +
  ylim(0, 5) +
  labs(x = "Particles/ml",
       y = "Length, mm")

#' we can use linear model to analyse length data
#' because we have some outliers, it's good idea to use robust model where
#' error is modeled by t distribution (not by normal distribution as default) 
#' 
library(MASS)
m1 <- rlm(len ~ conc, data = msml_filt)
summary(m1)
#' here the model is pretty clear that there is no effect on pod length
#' and slope is not different from zero and their average lenght is 4 mm

#' for comparison, regular linear model
m2 <- lm(len ~ conc, data = msml_filt)
summary(m2)

#' we have also count data 
#' here we use glm() function with poisson errors
m3 <- glm(n ~ conc, data = msml_filt, family = poisson)
#' seems that we have huge overdispersion: Residual deviance is much, much more bigger than degrees of freedom. This model is no good!
#' Residual deviance: 533.74  on 34  degrees of freedom
summary(m3)

#' there is even test for overdispersion
library(AER)
dispersiontest(m3, trafo = 1)

#' first we can try quasipoisson distribution
m4 <- glm(n ~ conc, data = msml_filt, family = quasipoisson)
summary(m4)
#' still no luck!

#' then, let's try to fit negative binomial model
m5 <- glm.nb(n ~ conc, data = msml_filt)
summary(m5)

#' predictions are identical...
library(modelr)
grid <- msml_filt %>% 
  gather_predictions(m3, m5) %>% 
  mutate(n = exp(pred))

ggplot(data = msml_filt, aes(x = log10(conc), y = n)) +
  geom_point() +
  stat_summary(fun.data = mean_se, color = "red") +
  geom_line(data = grid, aes(color = model)) +
  labs(x = "Particles/ml", y = "Number of pedigree")

#' so which model is better..
#' The Poisson and negative binomial (NB) model are nested: Poisson is a special case with theta = infinity. So a likelihood ratio test comparing the two models is testing the null hypothesis that "theta = infinity" against the alternative that "theta < infinity". Here the two models have the following log-likelihoods
logLik(m3) # poisson
logLik(m5) # neg. binom

#' neg. binomial fit seems way better, although predictions don't differ..
stat <- as.numeric(2 * (logLik(m5) - logLik(m1)))
pchisq(stat, df = 3 - 2)

library("lmtest")
lrtest(m3, m5)

#' one could also try to fit log-logistic model using poisson distribution with drc package
