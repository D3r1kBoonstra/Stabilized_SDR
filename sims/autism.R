# Packages ----------------------------------------------------------------
library("tidyverse");theme_set(theme_minimal());library("qqplotr")
library("hldr")
library("SCPME")
source("sim_wrappers.R")

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

set.seed(1)
dat[,2:11] <- dat[,2:11] + matrix(rnorm(nrow(dat[,2:11])*ncol(dat[,2:11]), sd = .00001), nrow = nrow(dat))

# {
#   x <- now()
#   omega <- SCPME_qda(dat[,-1], dat$class, nlam = 1000, cores = 7, K = 10, lam.max = 130)
#   lambdas <- vapply(seq_along(omega), function(i) omega[[i]]$Tuning[[2]], 
#                     FUN.VALUE = numeric(1))
#   now() - x
# }

## Saving Precision Estimates for eigen analysis
sim_autism_prec <- prec_est_sim(dat, 1000, lam = c(.86, 1.15))
# save(sim_autism_prec, file = "saved_sims/sim_autism_prec.RData")

## Qda shrink sims Repeated 10 fold CV
sim_autism_qda <- qda_shrink_sim(dat, 1000, lam = c(.86, 1.15))
# save(sim_autism_qda, file = "saved_sims/sim_autism_qda.RData")

### Median & se
load(file = "saved_sims/sim_autism_qda.RData")
sim_autism_qda |> 
  apply(2, function(x) paste0(round(median(x), 4)*100,"(", round(sd(x), 4)*100, ")"))


# HLDR Shrinkage ----------------------------------------------------------
## SY hldr shrink Repeated 10 fold CV
sim_autism_hldr <- hldr_sim(data = dat, lam = c(.73, .33), nsims = 1000)
# save(sim_autism_hldr, file = "saved_sims/sim_autism_hldr.RData")

## Median & SE
load(file = "saved_sims/sim_autism_hldr.RData")
lapply(sim_autism_hldr, function(x) {
  meds <- apply(x, 2 , median)
  ses <- apply(x, 2, sd)
  min_dim <- which.min(meds)
  paste0(round(meds[min_dim], 4)*100, "(", round(ses[min_dim], 5)*100, ")", "(", min_dim, ")")
}) |> 
  unlist()

sim_hldr_boxplot(sim_autism_hldr, dims = 1:9)

# HLDR Matrices -----------------------------------------------------------
## Saving HLDR Matrices for analysis
sim_autism_mats <- hldr_mats_sim(data = dat, nsims = 1000, lam =  c(.73, .33), dims = c(6, 5, 5, 4, 2)) 
# save(sim_autism_mats, file = "saved_sims/sim_autism_mats.RData")
