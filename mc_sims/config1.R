source("wrapper_fns.R")

# Configuration -----------------------------------------------------------
p <- 10
mu <- list(rep(0, p), rep(1, p))
Sigma <- list(diag(p), diag(p))


# Simulation --------------------------------------------------------------
config1 <- sim_wrapper(mu, Sigma, n = 5000, ni = c(p + 1, 2*p, 6*p), nsims = 1000, 
                       lambdas = c(1.45, .5))
# save(config1, file = "saved_sims/config1.RData")
load(file = "saved_sims/config1.RData")
config1 |> sim_summary()
config1 |> sim_boxplots()
