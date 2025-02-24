source("wrapper_fns.R")

# Configuration -----------------------------------------------------------
p <- 10
mu <- list(
  {
  mu1 <- c(-1.43, -0.66, -0.94, 0.31, -0.19, 0.89, 0.25, -0.34, 1.25, -1.60)
  mu1
  },
  mu1 + 1, mu1 + 2)

Sigma <- list(diag(p),
              diag(p) + 1,
              {
                S <- diag(p) + 1
                S[, 3] <- 0
                S[3, ] <- 0
                S[3, 3] <- 10
                S
              })

# Simulation --------------------------------------------------------------
config2 <- sim_wrapper(mu, Sigma, n = 5000, ni = c(p + 1, 2*p, 6*p), nsims = 1000, 
                       lambdas = c(.1, .06, .39), type = "qda")
# save(config2, file = "saved_sims/config2.RData")
load(file = "saved_sims/config2.RData")
config2 |> sim_summary()
config2 |> sim_boxplots()