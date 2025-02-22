# Packages ----------------------------------------------------------------
library("tidyverse");theme_set(theme_minimal());library("qqplotr")
library("hldr")
library("SCPME")
source("wrapper_fns.R")

#  Reading Data & Cleaning ------------------------------------------------
dat <- read_delim("datasets/divorce.csv", delim = ";") |> 
  janitor::clean_names() |> 
  relocate(class, .before = everything()) |>
  mutate(across(where(is.character), as_factor),
         across(-class, as.numeric),
         "class" = case_when(class  == "0" ~ 1,
                             class == "1" ~ 2) |>
           as_factor()
  )
set.seed(1)
dat[,-1] <- dat[,-1] + matrix(rnorm(nrow(dat[,-1])*ncol(dat[,-1]), sd = .00001), nrow = nrow(dat))

## Saving Precision Estimates for eigen analysis
sim_div_prec <- prec_est_sim(dat, 1000, lam = c(.2311, .031), type = "qda", gamma = c(1.5, 3))
# save(sim_div_prec, file = "saved_sims/sim_div_prec.RData")

## Qda shrink sims Repeated 10 fold CV
sim_div_qda <- qda_shrink_sim(dat, 1000, lam = c(.2311, .031), type = "qda", gamma = c(1.5, 3))
# save(sim_div_qda, file = "saved_sims/sim_div_qda.RData")

### Median & se
load(file = "saved_sims/sim_div_qda.RData")
sim_div_qda |> 
  apply(2, function(x) paste0(round(median(x), 4)*100,"(", round(sd(x), 4)*100, ")"))

# HLDR Shrinkage ----------------------------------------------------------

## SY hldr shrink Repeated 10 fold CV
sim_div_hldr <- hldr_sim(data = dat, lam = c(1.5, 1.5), nsims = 10, 
                         type = "qda", gamma = c(.55, 3.04))
# save(sim_div_hldr, file = "saved_sims/sim_div_hldr.RData")

load(file = "saved_sims/sim_div_hldr.RData")
lapply(sim_div_hldr, function(x) {
  meds <- apply(x, 2 , median)
  ses <- apply(x, 2, sd)
  min_dim <- which.min(meds)
  paste0(round(meds[min_dim], 4)*100, "(", round(ses[min_dim], 5)*100, ")", "(", min_dim, ")")
}) |> 
  unlist()
sim_hldr_boxplot(sim_div_hldr, dims = 2:9)

sim_div_comp <- comp_sim(dat, ends_u = 1, nsims = 10)
load(file = "saved_sims/sim_div_comp.RData")
apply(sim_div_comp, 2, 
      function(y) paste0(round(median(y), 4), "(", round(sd(y), 4), ")"))
# save(sim_div_comp, file = "saved_sims/sim_div_comp.RData")
(apply(sim_div_comp, 2, median)*100) |> round(2)

# HLDR Matrices -----------------------------------------------------------
## Saving HLDR Matrices for analysis
sim_div_mats <- hldr_mats_sim(data = dat, lam = c(1.5, 1.5), nsims = 1000, 
                              type = "qda", gamma = c(.55, 3.04), dims = c(2, 2 , 3, 1, 3)) 
# save(sim_div_mats, file = "saved_sims/sim_div_mats.RData")
