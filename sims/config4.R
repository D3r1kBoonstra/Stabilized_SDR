library("tidyverse");theme_set(theme_minimal())
library("sdr")
source("wrapper_fns.R")

# Configuration -----------------------------------------------------------
p <- 50
set.seed(1)
mu <- list(rep(0, p), runif(p, -.5, .5))
Sigma <- list(
  diag(p), diag(p) + 2
)

# Simulation --------------------------------------------------------------
config4 <- sim_wrapper(mu, Sigma, n = 5000, ni = c(p + 1, 2*p, 6*p), nsims = 1000,
                       lambdas = c(1, .01), type = "qda", gamma = c(1.4, 1.4))
# save(config4, file = "saved_sims/config4.RData")
load(file = "saved_sims/config4.RData")
config4 |> sim_summary()
config4 |> sim_boxplots()


