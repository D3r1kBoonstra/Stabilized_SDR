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

# Comparing Prec Est in SDRS ----------------------------------------------------------
sim_bc_sdrs <- sdrs_sim(data = dat, lam = c(3.5, .15), nsims = 1000, gamma = c(2.28, 4.9), 
                        type = "qda", standardize_xbar = TRUE)
# save(sim_bc_sdrs, file = "saved_sims/sim_bc_sdrs.RData")

## Median & SE
load(file = "saved_sims/sim_bc_sdrs.RData")
lapply(sim_bc_sdrs, function(x) {
  meds <- apply(x, 2 , median)
  ses <- apply(x, 2, sd)
  min_dim <- which.min(meds)
  paste0(round(meds[min_dim], 5)*100, "(", round(ses[min_dim], 5)*100, ")", "(", min_dim, ")")
}) |> 
  unlist()

## Boxplots
sim_sdrs_boxplot(sim_bc_sdrs, dims = 2:8)

# Competitors -------------------------------------------------------------
sim_bc_comp <- comp_sim(dat, ends_u = 1, nsims = 1000)
# save(sim_bc_comp, file = "saved_sims/sim_bc_comp.RData")
load(file = "saved_sims/sim_bc_comp.RData")
apply(sim_bc_comp, 2, 
      function(y) paste0(round(median(y), 4), "(", round(sd(y), 4), ")"))

