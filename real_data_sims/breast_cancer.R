# Packages ----------------------------------------------------------------
library("tidyverse");theme_set(theme_minimal());library("qqplotr")
library("hldr")
library("SCPME")
source("wrapper_fns.R")

#  Reading Data & Cleaning ------------------------------------------------
{
  data(BreastCancer, package = "mlbench")
  dat <- BreastCancer |> 
    select(-Id) |> 
    janitor::clean_names() |> 
    mutate(across(-class, as.numeric), 
           "class" = case_when(class == "benign" ~ 1, 
                               class == "malignant"~ 2)) |> 
    drop_na() |> 
    relocate(class, .before = everything()) 
  rm(BreastCancer)
}
## BC hldr benfits from dimselect

# Precision Estimation ----------------------------------------------------
## Saving Precision Estimates for eigen analysis
sim_bc_prec <- prec_est_sim(dat, 1000, lam = c(3.5, .15), gamma = c(2.28, 4.9), 
                            type = "qda", standardize_xbar = TRUE)
# save(sim_bc_prec, file = "saved_sims/sim_bc_prec.RData")
load(file = "saved_sims/sim_bc_prec.RData")
prec_analysis(sim_bc_prec)

scales::trans
# QDF Shrinkage -----------------------------------------------------------
## qda shrinkage repeated 10 fold CV
sim_bc_qda <- qda_shrink_sim(dat, 1000, lam = c(3.5, .15), gamma = c(2.28, 4.9), 
                             type = "qda", standardize_xbar = TRUE)
# save(sim_bc_qda, file = "saved_sims/sim_bc_qda.RData")

## Median & se
load(file = "saved_sims/sim_bc_qda.RData")
sim_bc_qda |> 
  apply(2, function(x) paste0(round(median(x), 4)*100,"(", round(sd(x), 5)*100, ")"))

# HLDR Shrinkage ----------------------------------------------------------
## hldr shrinkage repeated 10 fold CV
sim_bc_hldr <- hldr_sim(data = dat, lam = c(3.5, .15), nsims = 10, gamma = c(2.28, 4.9), 
                        type = "qda", standardize_xbar = TRUE)
# save(sim_bc_hldr, file = "saved_sims/sim_bc_hldr.RData")

## Median & se
load(file = "saved_sims/sim_bc_hldr.RData")
lapply(sim_bc_hldr, function(x) {
  meds <- apply(x, 2 , median)
  ses <- apply(x, 2, sd)
  min_dim <- which.min(meds)
  paste0(round(meds[min_dim], 5)*100, "(", round(ses[min_dim], 5)*100, ")", "(", min_dim, ")")
}) |> 
  unlist()

## Boxplots
sim_hldr_boxplot(sim_bc_hldr, dims = 2:8)
# u = 1
sim_bc_comp <- comp_sim(dat, ends_u = 1, nsims = 1000)

# save(sim_bc_comp, file = "saved_sims/sim_bc_comp.RData")
load(file = "saved_sims/sim_bc_comp.RData")
apply(sim_bc_comp, 2, 
      function(y) paste0(round(median(y), 4), "(", round(sd(y), 4), ")"))

# HLDR Matrices -----------------------------------------------------------
## Saving HLDR Matrices for analysis
sim_bc_mats <- hldr_mats_sim(dat, 1000, lam =  c(3.5, .15), dims = c(2, 3, 2, 2, 2), 
                             gamma = c(2.28, 4.9), 
                             type = "qda", standardize_xbar = TRUE) 
save(sim_bc_mats, file = "saved_sims/sim_bc_mats.RData")

