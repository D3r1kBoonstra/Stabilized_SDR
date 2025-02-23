# Packages ----------------------------------------------------------------
library("tidyverse");theme_set(theme_minimal())
library("sdr")
source("wrapper_fns.R")

#  Reading Data & Cleaning ------------------------------------------------
dat <- read_csv("datasets/ecoli.csv", col_names = FALSE) |> 
  janitor::clean_names() |> 
  rename("class" = last_col()) |> 
  relocate(class, .before = everything()) |> 
  filter(!class %in% c("imL", "imS", "omL")) |> 
  mutate("class" = as_factor(class) |> as.numeric() |> as_factor())

## adding error for singularity issues
set.seed(1)
dat$x3 <- dat$x3 + rnorm(length(dat$x3), sd = .00001)
dat$x4 <- dat$x4 + rnorm(length(dat$x4), sd = .00001)

# Comparing Prec Est in SDRS ----------------------------------------------------------
sim_ecoli_sdrs <- sdrs_sim(data = dat, lam = c(.05, .04, .04, .03, .04), nsims = 1000)
# save(sim_ecoli_sdrs, file = "saved_sims/sim_ecoli_sdrs.RData")

## Median & SE
load(file = "saved_sims/sim_ecoli_sdrs.RData")
lapply(sim_ecoli_sdrs, function(x) {
  meds <- apply(x, 2 , median)
  ses <- apply(x, 2, sd)
  min_dim <- which.min(meds)
  paste0(round(meds[min_dim], 4)*100, "(", round(ses[min_dim], 5)*100, ")", "(", min_dim, ")")
}) |> 
  unlist()

## Boxplots
sim_sdrs_boxplot(sim_ecoli_sdrs, dims = 4:6)

# Competitors -------------------------------------------------------------
sim_ecoli_comp <- comp_sim(data = dat, ends_u = 5, sir_d = 5, nsims = 1000)
# save(sim_ecoli_comp, file = "saved_sims/sim_ecoli_comp.RData")
load(file = "saved_sims/sim_ecoli_comp.RData")
apply(sim_ecoli_comp, 2, 
      function(y) paste0(round(median(y), 4), "(", round(sd(y), 4), ")"))


