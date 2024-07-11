# Packages ----------------------------------------------------------------
library("tidyverse");theme_set(theme_minimal());library("qqplotr")
library("hldr")
library("SCPME")
source("sim_wrappers.R")

# Reading Data & Cleaning -------------------------------------------------
dat <- read_csv("datasets/bands.data", col_names = FALSE, na = "?") |> 
  janitor::clean_names() |> 
  rename("class" = last_col()) |>
  select(-c(1:4)) |> 
  relocate(class, .before = everything()) |>
  mutate(across(where(is.character), as_factor),
         across(-class, as.numeric), 
         "class" = case_when(class == "band" ~ 1, 
                             class == "noband" ~ 2) |> 
           as_factor()
  ) |> 
  drop_na()

set.seed(1)
dat[,-1] <- dat[,-1] + matrix(rnorm(nrow(dat[,-1])*ncol(dat[,-1]), sd = .00001), nrow = nrow(dat))
dat <- dat |> 
  mutate(across(where(is.numeric), function(x)(x - mean(x))/sd(x)))
## Getting SCPME tuning parameters
{
  x <- now()
  omega <- SCPME_qda(dat[,-1], dat$class, nlam = 1000, cores = 7, K = 10, lam.max = 120)
  lambdas <- vapply(seq_along(omega), function(i) omega[[i]]$Tuning[[2]], 
                    FUN.VALUE = numeric(1))
  now() - x
  }

## Saving Precision Estimates for eigen analysis
sim_bands_prec <- prec_est_sim(dat, 250, lam = c(.105, .106))
# save(sim_bands_prec, file = "saved_sims/sim_bands_prec.RData")

## Qda shrink sims Repeated 10 fold CV
sim_bands_qda <- qda_shrink_sim(dat, 20, lam = c(.10, .109), gamma = c(2, 2))
# save(sim_bands_qda, file = "saved_sims/sim_bands_qda.RData")

### Median & se
load(file = "saved_sims/sim_bands_qda.RData")
sim_bands_qda |> 
  apply(2, function(x) paste0(round(median(x, na.rm = TRUE), 4)*100,"(", round(sd(x), 4)*100, ")"))

# HLDR Matrices -----------------------------------------------------------
## Saving HLDR Matrices for analysis
sim_bands_mats <- hldr_mats_sim(dat, 250, lam =  c(.12, .12)) 
# save(sim_bands_mats, file = "saved_sims/sim_bands_mats.RData")

# HLDR Shrinkage ----------------------------------------------------------

## SY hldr shrink Repeated 10 fold CV
sim_bands_hldr <- hldr_sim(data = dat, lam = c(.12, .12), nsims = 250)
# save(sim_bands_hldr, file = "saved_sims/sim_bands_hldr.RData")

load(file = "saved_sims/sim_bands_hldr.RData")
lapply(sim_bands_hldr, function(x) {
  meds <- apply(x, 2 , median)
  ses <- apply(x, 2, sd)
  min_dim <- which.min(meds)
  paste0(round(meds[min_dim], 4)*100, "(", round(ses[min_dim], 5)*100, ")", "(", min_dim, ")")
}) |> 
  unlist()

sim_hldr_boxplot(sim_bands_hldr, dims = 10:32)
