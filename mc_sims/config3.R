library("tidyverse");theme_set(theme_minimal())
library("sdr")
source("wrapper_fns.R")


# Configuration -----------------------------------------------------------
p <- 10
mu <- list(rep(0, p), rep(5, p), rep(10, p))
Sigma <- list({
  S <- 
    matrix(c(10, 4, 5, 4, 3, 4, 4, 5, 4, 3, 
             4, 10, 5, 2, 4, 3, 3, 5, 4, 3, 
             5, 5, 10, 5, 5, 3, 4, 4, 4, 4, 
             4, 2, 5, 10, 3, 4, 2, 3, 4, 3, 
             3, 4, 5, 3, 12, 3, 4, 5, 3, 3, 
             4, 3, 3, 4, 3, 9, 3, 4, 4, 4, 
             4, 3, 4, 2, 4, 3, 14, 2, 2, 2, 
             5, 5, 4, 3, 5, 4, 2, 12, 1, -.5, 
             4, 4, 4, 4, 3, 4, 2, 1, 14, -1, 
             3, 3, 4, 3, 3, 4, 2, -.5, -1, 11), 
           nrow = 10, ncol = 10, byrow = TRUE)
}, 
S, 
diag(p))


# Simulation --------------------------------------------------------------
config3 <- sim_wrapper(mu, Sigma, n = 5000, ni = c(p + 1, 2*p, 6*p), nsims = 1000, 
                       lambdas = c(.03, .03, .13))
# save(config3, file = "saved_sims/config3.RData")
load(file = "saved_sims/config3.RData")
config3 |> sim_summary()
config3 |> sim_boxplots()
