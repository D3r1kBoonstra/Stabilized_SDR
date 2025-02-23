# Packages ----------------------------------------------------------------
library("tidyverse");theme_set(theme_minimal())
library("sdr")
source("wrapper_fns.R")

#  Reading Data & Cleaning ------------------------------------------------
dat <- read_csv("datasets/Autism-Adolescent-Data.arff", 
                na = "?", col_names = FALSE) |> 
  janitor::clean_names() |> 
  rename("class" = last_col()) |>
  select(-c(x17, x19)) |> 
  relocate(class, .before = everything()) |>
  mutate(across(where(is.character), as_factor),
         across(-class, as.numeric), 
         "class" = case_when(class == "NO" ~ 1, 
                             class == "YES" ~ 2) |> 
           as_factor()
  ) |> 
  drop_na()

## adding error for singularity issues
set.seed(1)
dat[,2:11] <- dat[,2:11] + matrix(rnorm(nrow(dat[,2:11])*ncol(dat[,2:11]), sd = .00001), nrow = nrow(dat))

# Comparing Prec Est in SDRS ----------------------------------------------------------
sim_autism_sdrs <- sdrs_sim(data = dat, lambdas = c(.73, .33), nsims = 1000)
# save(sim_autism_sdrs, file = "saved_sims/sim_autism_sdrs.RData")

## Median & SE
load(file = "saved_sims/sim_autism_sdrs.RData")
lapply(sim_autism_sdrs, function(x) {
  meds <- apply(x, 2 , median)
  ses <- apply(x, 2, sd)
  min_dim <- which.min(meds)
  paste0(round(meds[min_dim], 4)*100, "(", round(ses[min_dim], 5)*100, ")", "(", min_dim, ")")
}) |> 
  unlist()

## Boxplots
sim_sdrs_boxplot(sim_autism_sdrs, dims = 1:9)


# Competitors -------------------------------------------------------------
sim_autism_comp <- comp_sim(dat, ends_u = 2, nsims = 1000)
# save(sim_autism_comp, file = "saved_sims/sim_autism_comp.RData")
load(file = "saved_sims/sim_autism_comp.RData")
apply(sim_autism_comp, 2, function(y) paste0(round(median(y), 4), "(", round(sd(y), 4), ")"))
