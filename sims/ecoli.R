# Packages ----------------------------------------------------------------
library("tidyverse");theme_set(theme_minimal());library("qqplotr")
library("hldr")
library("SCPME")
source("sim_wrappers.R")

#  Reading Data & Cleaning ------------------------------------------------
dat <- read_csv("datasets/ecoli.csv", col_names = FALSE) |> 
  janitor::clean_names() |> 
  rename("class" = last_col()) |> 
  relocate(class, .before = everything()) |> 
  filter(!class %in% c("imL", "imS", "omL")) |> 
  mutate("class" = as_factor(class))

set.seed(1)
dat$x3 <- dat$x3 + rnorm(length(dat$x3), sd = .00001)
dat$x4 <- dat$x4 + rnorm(length(dat$x4), sd = .00001)

## Saving Precision Estimates for eigen analysis
sim_ecoli_prec <- prec_est_sim(dat, 1000, lam = c(.05, .04, .04, .03, .04))
# save(sim_ecoli_prec, file = "saved_sims/sim_ecoli_prec.RData")

sim_ecoli_qda <- qda_shrink_sim(dat, 1000, lam = c(.05, .04, .04, .03, .04))
# save(sim_ecoli_qda, file = "saved_sims/sim_ecoli_qda.RData")

## Median & SE
load(file = "saved_sims/sim_ecoli_qda.RData")
sim_ecoli_qda |> 
  apply(2, function(x) paste0(round(median(x), 4)*100,"(", round(sd(x), 5)*100, ")"))

# HLDR Matrices -----------------------------------------------------------
## Saving HLDR Matrices for analysis
sim_ecoli_mats <- hldr_mats_sim(data = dat, nsims = 1000, 
                                lam =  c(.05, .04, .04, .03, .04), dims = c(6, 6, 5, 5, 4)) 
# save(sim_ecoli_mats, file = "saved_sims/sim_ecoli_mats.RData")

# HLDR Shrinkage ----------------------------------------------------------
sim_ecoli_hldr <- hldr_sim(data = dat, lam = c(.05, .04, .04, .03, .04), nsims = 1000)
# save(sim_ecoli_hldr, file = "saved_sims/sim_ecoli_hldr.RData")

load(file = "saved_sims/sim_ecoli_hldr.RData")
lapply(sim_ecoli_hldr, function(x) {
  meds <- apply(x, 2 , median)
  ses <- apply(x, 2, sd)
  min_dim <- which.min(meds)
  paste0(round(meds[min_dim], 4)*100, "(", round(ses[min_dim], 5)*100, ")", "(", min_dim, ")")
}) |> 
  unlist()

sim_hldr_boxplot(sim_ecoli_hldr, dims = 4:6)


