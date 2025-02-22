# Packages ----------------------------------------------------------------
library("tidyverse");theme_set(theme_minimal());library("qqplotr")
library("hldr")
library("SCPME")
source("wrapper_fns.R")

#  Reading Data & Cleaning ------------------------------------------------
dat <- read_csv("datasets/wheat-seeds.csv", 
                na = "?", col_names = FALSE) |> 
  janitor::clean_names() |> 
  rename("class" = last_col()) |>
  relocate(class, .before = everything()) |>
  mutate("class" = as_factor(class)) |> 
  drop_na() 


# Precision Estimation ----------------------------------------------------
## Saving Precision Estimates for eigen analysis
sim_wheat_prec <- prec_est_sim(dat, 1000, lam = c(.00021, .00011, .00021))
save(sim_wheat_prec, file = "saved_sims/sim_wheat_prec.RData")


# QDF Shrinkage -----------------------------------------------------------
sim_wheat_qda <- qda_shrink_sim(dat, 1000, lam = c(.00021, .00011, .00021))
# save(sim_wheat_qda, file = "saved_sims/sim_wheat_qda.RData")

## Median & SE
load(file = "saved_sims/sim_wheat_qda.RData")
sim_wheat_qda |> 
  apply(2, function(x) paste0(round(median(x), 4)*100,"(", round(sd(x), 4)*100, ")"))

# HLDR Shrinkage ----------------------------------------------------------
## hldr shrinkage repeated 10 fold CV
sim_wheat_hldr <- hldr_sim(data = dat, lam = c(0.00159, 0.00205, 0.00154), nsims = 1000)
# save(sim_wheat_hldr, file = "saved_sims/sim_wheat_hldr.RData")

load(file = "saved_sims/sim_wheat_hldr.RData")
lapply(sim_wheat_hldr, function(x) {
  meds <- apply(x, 2 , median)
  ses <- apply(x, 2, sd)
  min_dim <- which.min(meds)
  paste0(round(meds[min_dim], 4)*100, "(", round(ses[min_dim], 5)*100, ")", "(", min_dim, ")")
}) |> 
  unlist()
## Got .0190 once but cant get again

sim_hldr_boxplot(sim_wheat_hldr)

sim_wheat_comp <- comp_sim(data = dat, ends_u = 4, sir_d = 2, nsims = 1000)
load(file = "saved_sims/sim_wheat_comp.RData")
apply(sim_wheat_comp, 2, 
      function(y) paste0(round(median(y), 4), "(", round(sd(y), 4), ")"))
# save(sim_wheat_comp, file = "saved_sims/sim_wheat_comp.RData")
(apply(sim_wheat_comp, 2, median)*100) |> round(2)

# HLDR Matrices -----------------------------------------------------------
## Saving HLDR Matrices for analysis
sim_wheat_mats <- hldr_mats_sim(dat, 1000, lam =  c(0.00159, 0.00205, 0.00154), dims = c(rep(6, 4), 5)) 
# save(sim_wheat_mats, file = "saved_sims/sim_wheat_mats.RData")




