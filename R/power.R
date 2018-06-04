
library(tidyverse)
library(pwr)

# Our null hypothesis is that the coin is fair and lands heads 50% of the time (π=0.50)
# Our alternative hypothesis is that the coin is loaded to land heads more then 50% of the time (π>0.50).
?pwr.p.test
?ES.h
(p.out <- pwr.p.test(h = ES.h(p1 = 0.75, p2 = 0.50), 
           sig.level = 0.05, 
           power = 0.80, 
           alternative = "greater"))

plot(p.out)


# base R power
?power.anova.test()
groupmeans <- c(120, 130, 140, 150)
power.anova.test(groups = length(groupmeans),
                 between.var = var(groupmeans),
                 within.var = 500, 
                 power = .8)

# anova 2x2 power

sqrt(16 / 6)
sqrt(16) / sqrt(6)
