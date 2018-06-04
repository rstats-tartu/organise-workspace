
## ----------------------- ##
## Computing and Graphing  ##
## Statistical Power       ##
## ----------------------- ##

## Parameters of the Problem ##

mu0 <- 6     ## Null-hypothesized value
alpha <- 0.05  ## Significance Level
s.sq <- 16    ## Variance from pilot study
n1 <- 10    ## Sample Size
alt <- 6.5

se <- sqrt(s.sq / n1 - 1)

## -------------------------- ##
## What is Statistical Power? ##
## -------------------------- ##

cuts = c(1 - alpha / 2, alpha / 2)
crits = qnorm(cuts, mu0, sqrt(s.sq / n1))    ## Critical Values Based on Null
shadenorm(mu = mu0, sig = se, outside = crits) ## My own function
shadenorm(mu = 6.5, sig = se, lines = TRUE, outside = crits, col = "blue")

## ------------------------------ ##
## Write a Function to Compute it ##
## ------------------------------ ##

#' @param theta alternative hypothesis
#' @param mu null hypothesis
#' @param var variance
#' @param n population or sample size
#' @param alpha significance level 
power <- function(theta, mu, var, n, alpha = 0.05) {
  crit_lo <- qnorm(alpha / 2, mu, sqrt(var / n))    
  crit_up <- qnorm(1 - alpha / 2, mu, sqrt(var / n))  
  pr_up <- pnorm(crit_up, theta, sd = sqrt(var / n), lower.tail = FALSE)
  pr_lo <- pnorm(crit_lo, theta, sd = sqrt(var / n))
  pr_lo + pr_up
}

power(theta = 9, mu = 6, var = 16, n = 14)
theta <- seq(1, 12, by = 0.01)
pow <- power(theta, mu0, s.sq, n1)
plot(theta, pow, type = "l", ylim = c(0, 1), main = "My Power Plot")

