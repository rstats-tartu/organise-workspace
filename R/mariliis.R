
#' check and install required packages
install.packages("pacman")
pacman::p_load(tidyverse, readxl, car, drc, sandwich, lmtest, multcomp)

#' load required libraries
library(tidyverse)
library(readxl)

#' import data file from excel file
#' get sheet names
(sheets <- excel_sheets("docs/Statistika_Küsimused.xlsx"))
#' import region of interest from sheet
mariliis <- read_excel("docs/Statistika_Küsimused.xlsx", 
                       sheet = "Mariliis",
                       range = c("A10:E34"))
#' fix imported table - 
#' remove bigger than mark, 
#' replace comma with period, 
#' convert values to numbers
#' fix column names
mariliis_fix <- mutate_at(mariliis, vars(starts_with("EC50")), str_remove, pattern = ">") %>% 
  mutate_at(vars(starts_with("EC50")), str_replace, pattern = ",", replacement = "\\.") %>% 
  mutate_at(vars(starts_with("EC50")), parse_number) %>% 
  dplyr::select(c_ahela_pikkus = `C-ahela pikkus (C2-C16)`,
         katioon = `Katioon: Py, Imid või Chol`,
         ec50_bak1 = `EC50 (µM), Bakter1`,
         ec50_bak2 = `EC50 (µM), Bakter2`)
#'NB! dplyr::select() becomes masked by MASS::select(), therefore we need to use here dplyr::select(). Otherwise you can use just select().


#' how does ec50 values look against normal distribution
library(car)
qqPlot(c(mariliis_fix$ec50_bak1, mariliis_fix$ec50_bak2))
#' values don't seem to follow normal distribution

#' to test two bacteria using Wilcoxon rank sum test
#' wilcoxon tests medians 
wilcox.test(mariliis_fix$ec50_bak2, mariliis_fix$ec50_bak1, paired = TRUE)

#' let's take step back
(mariliis_tidy <- gather(mariliis_fix, 
                        key = species, 
                        value = ec50, 
                        -c_ahela_pikkus, 
                        -katioon) %>% 
  mutate(species = str_remove(species, "ec50_")))

#' calculate summary statistics for dataset
group_by(mariliis_tidy, species) %>%
  summarise(count = n(),
            median = median(ec50, na.rm = TRUE),
            IQR = IQR(ec50, na.rm = TRUE))

#' make a boxplot of data
ggplot(data = mariliis_tidy, mapping = aes(x = species, y = ec50)) +
  geom_boxplot() +
  geom_jitter(width = 0.1)

#' formula version of wilcoxon test
wilcox.test(ec50 ~ species, 
            data = mariliis_tidy, 
            paired = TRUE)
#' note that V = 0 now. The value V corresponds to the sum of ranks assigned 
#' to the differences with positive sign. 
#' We can try to filter out values that are equal between two bacteria:
#' filter(mariliis_tidy, ec50 != 2000)

#' Next, we are interested what chemicals affect toxicity most
(m1 <- lm(ec50 ~ katioon + c_ahela_pikkus, data = mariliis_tidy))
#' linear model indicates that compounds with c-chain length over 6 have higher toxicity
summary(m1)
#' anova analysis suggest that toxicity is indeed most influenced by c-chain length
(a1 <- aov(m1))
summary(a1)

#' calculating ec50 by using drc package
#' import inhibition data
(inh <- read_excel("docs/Statistika_Küsimused.xlsx", 
                       sheet = "Mariliis",
                       range = c("A44:E52")))

#' fix table, here we need to rename columns
inh_renamed <- rename(inh, conc = `Nominal concentration, mg/l`)
inh_tidy <- gather(inh_renamed, starts_with("EXPERIMENT"), key = "exp", value = "inhib") %>% 
  mutate(inhib = as.numeric(inhib))

ggplot(data = inh_tidy) +
  geom_point(mapping = aes(x = conc, y = inhib))

library(drc)
#' f(x) = c + \frac{d-c}{1+\exp(b(\log(x)-\log(e)))}
#' b is the steepness of the dose-response curve
#' c, d the lower and upper asymptotes/limits of the response
#' e the effective dose ED50
#' here we use four parameter model and constrain lower bound to 0,
#' we also constrain upper bound to 100
#' same can be achieved with LL.3(fixed = c(NA, 100, NA))
mod <- drm(inhib ~ conc, 
           data = inh_tidy, 
           fct = LL.4(fixed = c(NA, 0, 100, NA), names = c("steepness", 
                                                           "lower",
                                                           "upper",
                                                           "EC50")))
#' our coefs are pretty close to ones shown in origignal excele table
summary(mod)


#' let's plot fit for sanity check
plot(mod, type = "all",
     xlab = "Nominal conc, mg/l", xlim = c(0, 100),
     ylab = "Inhibition %")

#' we use the R packages lmtest and sandwich to obtain robust standard errors
library(sandwich)
library(lmtest)
coeftest(mod, vcov = sandwich)

#' By using simultaneous inference with function glht() in the R package
#' multcomp we can test wether coefficients differ from 0
library(multcomp)
summary(glht(mod))

#' ED() function can be used to estimate effective doses ED5, ED10, and ED50 #' is accomplished using 
ED(mod, c(5, 10, 50, 90), interval = "delta")

#' let's make high quality plot
# new dose levels as support for the line
newdata <- expand.grid(conc = exp(seq(log(0.5), log(100), length = 100)))

#' predictions and confidence intervals
pm <- predict(mod, newdata = newdata, interval = "confidence")
mod_pred <- as_data_frame(cbind(newdata, pm))

#' compose plot
ggplot(data = mod_pred) +
  geom_line(mapping = aes(x = conc, y = Prediction, group = 1)) +
  geom_ribbon(mapping = aes(x = conc, ymin = Lower, ymax = Upper), alpha = 0.2) +
  geom_point(data = inh_tidy, mapping = aes(x = conc, y = inhib)) +
  coord_trans(x = "log") +
  labs(x = "Nominal concentration, mg/l",
       y = "Luminescence inhibition, %",
       caption = "Line is four parameter log-logistic model fit\nwith lower and upper bounds fixed to 0 and 100.\nGray ribbon is 95% confidence interval for model fit.")
