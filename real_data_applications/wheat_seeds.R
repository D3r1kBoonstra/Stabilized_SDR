source("wrapper_fns.R")

#  Reading Data & Cleaning ------------------------------------------------
dat <- read_csv("datasets/wheat-seeds.csv", 
                na = "?", col_names = FALSE) |> 
  janitor::clean_names() |> 
  rename("class" = last_col()) |>
  relocate(class, .before = everything()) |>
  mutate("class" = as_factor(class)) |> 
  drop_na() 

# Comparing Prec Est in SDRS ----------------------------------------------------------
sim_wheat_sdrs <- sdrs_sim(data = dat, lam = c(0.00159, 0.00205, 0.00154), nsims = 1000)
# save(sim_wheat_sdrs, file = "saved_sims/sim_wheat_sdrs.RData")

## Median & SE
load(file = "saved_sims/sim_wheat_sdrs.RData")
lapply(sim_wheat_sdrs, function(x) {
  meds <- apply(x, 2 , median)
  ses <- apply(x, 2, sd)
  min_dim <- which.min(meds)
  paste0(round(meds[min_dim], 4)*100, "(", round(ses[min_dim], 5)*100, ")", "(", min_dim, ")")
}) |> 
  unlist()

## Boxplots
sim_sdrs_boxplot(sim_wheat_sdrs)

# Competitors -------------------------------------------------------------
sim_wheat_comp <- comp_sim(data = dat, ends_u = 4, sir_d = 2, nsims = 1000)
# save(sim_wheat_comp, file = "saved_sims/sim_wheat_comp.RData")
load(file = "saved_sims/sim_wheat_comp.RData")
apply(sim_wheat_comp, 2, 
      function(y) paste0(round(median(y), 4), "(", round(sd(y), 4), ")"))



