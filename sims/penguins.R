# Packages ----------------------------------------------------------------
library("tidyverse");theme_set(theme_minimal());library("qqplotr")
library("hldr")
library("SCPME")
source("sim_wrappers.R")

#  Reading Data & Cleaning ------------------------------------------------
dat <- palmerpenguins::penguins |> 
  janitor::clean_names() |> 
  select(-year) |> 
  rename("class" = species) |> 
  mutate(across(where(is.character), as_factor), 
         across(-class, as.numeric)) |> 
  drop_na()

set.seed(1)
dat$island <- dat$island + rnorm(length(dat$island), sd = .00001)

# {
#   x <- now()
#   omega <- SCPME_qda(dat[,-1], dat$class, nlam = 1000, cores = 7, K = 10, lam.max = 200)
#   lambdas <- vapply(seq_along(omega), function(i) omega[[i]]$Tuning[[2]], 
#                     FUN.VALUE = numeric(1))
#   now() - x
#   }

## Saving Precision Estimates for eigen analysis
sim_peng_prec <- prec_est_sim(dat, 1000, lam = c(.4, .1, .05), 
                              type = "qda", standardize_xbar = TRUE)
# save(sim_peng_prec, file = "saved_sims/sim_peng_prec.RData")

sim_peng_qda <- qda_shrink_sim(dat, 1000, lam = c(.4, .1, .05), 
                               type = "qda", standardize_xbar = TRUE)
# save(sim_peng_qda, file = "saved_sims/sim_peng_qda.RData")

### Median & se
load(file = "saved_sims/sim_peng_qda.RData")
sim_peng_qda |> 
  apply(2, function(x) paste0(round(median(x), 4)*100,"(", round(sd(x), 5)*100, ")"))


# HLDR Matrices -----------------------------------------------------------
## Saving HLDR Matrices for analysis
sim_peng_mats <- hldr_mats_sim(data = dat, lambdas = c(.4, .01, .008), nsims = 1000, 
                               type = "qda", standardize_xbar = TRUE, 
                               dims = c(5, 5, 5, 4, 5)) 
# save(sim_peng_mats, file = "saved_sims/sim_peng_mats.RData")

# HLDR Shrinkage ----------------------------------------------------------
sim_peng_hldr <- hldr_sim(data = dat, lam = c(.4, .01, .008), nsims = 1000, 
                          type = "qda", standardize_xbar = TRUE)
# save(sim_peng_hldr, file = "saved_sims/sim_peng_hldr.RData")

load(file = "saved_sims/sim_peng_hldr.RData")
lapply(sim_peng_hldr, function(x) {
  meds <- apply(x, 2 , median)
  ses <- apply(x, 2, sd)
  min_dim <- which.min(meds)
  paste0(round(meds[min_dim], 5)*100, "(", round(ses[min_dim], 5)*100, ")", "(", min_dim, ")")
}) |> 
  unlist()

sim_hldr_boxplot(sim_peng_hldr, dims = 3:5)+
  ylim(0, .15)


