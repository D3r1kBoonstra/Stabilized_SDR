# Packages ----------------------------------------------------------------
library("tidyverse");theme_set(theme_minimal());library("qqplotr")
library("hldr")
library("SCPME")
source("sim_wrappers.R")

#  Reading Data & Cleaning ------------------------------------------------
dat <- imputeR::spect |> 
  janitor::clean_names() |> 
  rename("class" = x1) |> 
  mutate("class" = as_factor(class))

set.seed(1)
dat[,-1] <- dat[,-1] + matrix(rnorm(nrow(dat[,-1])*ncol(dat[,-1]), sd = .00001), nrow = nrow(dat))

# {
#   x <- now()
#   omega <- SCPME_qda(dat[,-1], dat$class, nlam = 1000, cores = 7, K = 10)
#   lambdas <- vapply(seq_along(omega), function(i) omega[[i]]$Tuning[[2]], 
#                     FUN.VALUE = numeric(1))
#   now() - x
#   }

## Saving Precision Estimates for eigen analysis
sim_spect_prec <- prec_est_sim(dat, 1000, lam = c(.098, .27))
# save(sim_spect_prec, file = "saved_sims/sim_spect_prec.RData")

sim_spect_qda <- qda_shrink_sim(dat, 1000, lam = c(.098, .27))
# save(sim_spect_qda, file = "saved_sims/sim_spect_qda.RData")


### Median & se
load(file = "saved_sims/sim_spect_qda.RData")
sim_spect_qda |> 
  apply(2, function(x) paste0(round(median(x), 4)*100,"(", round(sd(x), 5)*100, ")"))


# HLDR Matrices -----------------------------------------------------------
## Saving HLDR Matrices for analysis
sim_spect_mats <- hldr_mats_sim(dat, 250, lam =  c(.045, .263), dims = c(21, 18, 2, 1, 1)) 
# save(sim_spect_mats, file = "saved_sims/sim_spect_mats.RData")

# HLDR Shrinkage ----------------------------------------------------------
## SY hldr shrink Repeated 10 fold CV
sim_spect_hldr <- hldr_sim(data = dat, lam = c(.045, .263), nsims = 1000)
# save(sim_spect_hldr, file = "saved_sims/sim_spect_hldr.RData")

load(file = "saved_sims/sim_spect_hldr.RData")
lapply(sim_spect_hldr, function(x) {
  meds <- apply(x, 2 , median)
  ses <- apply(x, 2, sd)
  min_dim <- which.min(meds)
  paste0(round(meds[min_dim], 4)*100, "(", round(ses[min_dim], 5)*100, ")", "(", min_dim, ")")
}) |> 
  unlist()

sim_hldr_boxplot(sim_spect_hldr)

