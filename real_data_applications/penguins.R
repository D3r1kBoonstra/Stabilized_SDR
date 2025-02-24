source("wrapper_fns.R")

#  Reading Data & Cleaning ------------------------------------------------
dat <- palmerpenguins::penguins |> 
  janitor::clean_names() |> 
  select(-year) |> 
  rename("class" = species) |> 
  mutate(across(where(is.character), as_factor), 
         across(-class, as.numeric), 
         "class" = as_factor(as.numeric(class))) |> 
  drop_na()

## adding error for singularity issues
set.seed(1)
dat$island <- dat$island + rnorm(length(dat$island), sd = .00001)

# Comparing Prec Est in SDRS ----------------------------------------------------------
sim_peng_sdrs <- sdrs_sim(data = dat, lam = c(.4, .01, .008), nsims = 1000, 
                          type = "qda", standardize_xbar = TRUE)
# save(sim_peng_sdrs, file = "saved_sims/sim_peng_sdrs.RData")

## Median & SE
load(file = "saved_sims/sim_peng_sdrs.RData")
lapply(sim_peng_sdrs, function(x) {
  meds <- apply(x, 2 , median)
  ses <- apply(x, 2, sd)
  min_dim <- which.min(meds)
  paste0(round(meds[min_dim], 5)*100, "(", round(ses[min_dim], 5)*100, ")", "(", min_dim, ")")
}) |> 
  unlist()

## Boxplots
sim_sdrs_boxplot(sim_peng_sdrs, dims = 3:5)


# Competitors -------------------------------------------------------------
sim_peng_comp <- comp_sim(dat, ends_u = 3,sir_d = 3,nsims = 1000)
# save(sim_peng_comp, file = "saved_sims/sim_peng_comp.RData")
load(file = "saved_sims/sim_peng_comp.RData")
apply(sim_peng_comp, 2, 
      function(y) paste0(round(median(y), 4), "(", round(sd(y), 4), ")"))

