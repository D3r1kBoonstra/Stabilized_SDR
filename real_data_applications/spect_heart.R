# Packages ----------------------------------------------------------------
library("tidyverse");theme_set(theme_minimal())
library("sdr")
source("wrapper_fns.R")

#  Reading Data & Cleaning ------------------------------------------------
dat <- imputeR::spect |> 
  janitor::clean_names() |> 
  rename("class" = x1) |> 
  mutate("class" = as_factor(class))

## adding error for singularity issues
set.seed(1)
dat[,-1] <- dat[,-1] + matrix(rnorm(nrow(dat[,-1])*ncol(dat[,-1]), sd = .00001), nrow = nrow(dat))

# Comparing Prec Est in SDRS ----------------------------------------------------------
sim_spect_sdrs <- sdrs_sim(data = dat, lam = c(.045, .263), nsims = 1000)
# save(sim_spect_sdrs, file = "saved_sims/sim_spect_sdrs.RData")

## Median & SE
load(file = "saved_sims/sim_spect_sdrs.RData")
lapply(sim_spect_sdrs, function(x) {
  meds <- apply(x, 2 , median)
  ses <- apply(x, 2, sd)
  min_dim <- which.min(meds)
  paste0(round(meds[min_dim], 4)*100, "(", round(ses[min_dim], 5)*100, ")", "(", min_dim, ")")
}) |> 
  unlist()

## Boxplots
sim_sdrs_boxplot(sim_spect_sdrs)


# Competitors -------------------------------------------------------------
sim_spect_comp <- comp_sim(dat, ends_u = 1, nsims = 1000)
# save(sim_spect_comp, file = "saved_sims/sim_spect_comp.RData")
load(file = "saved_sims/sim_spect_comp.RData")
apply(sim_spect_comp, 2, 
      function(y) paste0(round(median(y), 4), "(", round(sd(y), 4), ")"))

