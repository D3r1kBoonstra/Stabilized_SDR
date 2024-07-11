# Packages ----------------------------------------------------------------
library("tidyverse");theme_set(theme_minimal());library("qqplotr")
library("hldr")
library("SCPME")
source("sim_wrappers.R")

#  Reading Data & Cleaning ------------------------------------------------
dat <- read_csv("datasets/ionosphere.csv", col_names = FALSE, na = "?") |> 
  janitor::clean_names() |> 
  rename("class" = last_col()) |>
  select(-c(1:2)) |> 
  relocate(class, .before = everything()) |>
  mutate(across(where(is.character), as_factor),
         across(-class, as.numeric), 
         "class" = case_when(class == "g" ~ 1, 
                             class == "b" ~ 2) |> 
           as_factor()
  ) |> 
  drop_na() 

# Precision Estimation ----------------------------------------------------
## Saving Precision Estimates for eigen analysis
sim_ion_prec <- prec_est_sim(dat, 1000, lam = c(.0005, .0001), type = "qda",
                            gamma = c(1, 1.1))
# save(sim_ion_prec, file = "saved_sims/sim_ion_prec.RData")


# QDF Shrinkage -----------------------------------------------------------
## qda shrinkage repeated 10 fold CV
sim_ion_qda <- qda_shrink_sim(dat, 1000, lam = c(.0005, .0001), type = "qda",
                             gamma = c(1, 1.1))
# save(sim_ion_qda, file = "saved_sims/sim_ion_qda.RData")

## Median & se
load(file = "saved_sims/sim_ion_qda.RData")
sim_ion_qda |> 
  apply(2, function(x) paste0(round(median(x), 4)*100,"(", round(sd(x), 5)*100, ")"))

# HLDR Shrinkage ----------------------------------------------------------
## hldr shrinkage repeated 10 fold CV
sim_ion_hldr <- hldr_sim(data = dat, lam = c(.0005, .0005), nsims = 1000, 
                        gamma = c(1.04, .43), type = "qda", 
                        standardize_xbar = FALSE)
# save(sim_ion_hldr, file = "saved_sims/sim_ion_hldr.RData")

## Median & se
load(file = "saved_sims/sim_ion_hldr.RData")
lapply(sim_ion_hldr, function(x) {
  meds <- apply(x, 2 , median)
  ses <- apply(x, 2, sd)
  min_dim <- which.min(meds)
  paste0(round(meds[min_dim], 5)*100, "(", round(ses[min_dim], 5)*100, ")", "(", min_dim, ")")
}) |> 
  unlist()

## Boxplots
sim_hldr_boxplot(sim_ion_hldr, dims = 2:11)

# HLDR Matrices -----------------------------------------------------------
## Saving HLDR Matrices for analysis
sim_ion_mats <- hldr_mats_sim(data = dat, lam = c(.0005, .0005), nsims = 1000, 
                             gamma = c(1.04, .43), type = "qda", 
                             standardize_xbar = FALSE, dims = c(8, 8, 8, 8, 9)) 
# save(sim_ion_mats, file = "saved_sims/sim_ion_mats.RData")
