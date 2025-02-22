# Packages ----------------------------------------------------------------
library("tidyverse");theme_set(theme_minimal());library("qqplotr")
library("hldr")
library("SCPME")
source("wrapper_fns.R")

#  Reading Data & Cleaning ------------------------------------------------
dat <- read_csv("datasets/Dry_Bean_Dataset.arff", col_names = FALSE) |> 
  janitor::clean_names() |>
  rename("class" = last_col()) |> 
  relocate(class, .before = everything()) |>
  mutate(across(where(is.numeric), function(x)(x - mean(x))/sd(x)), 
         "class" = as_factor(class) |> as.numeric() |> as_factor())

## Saving Precision Estimates for eigen analysis
sim_beans_prec <- prec_est_sim(dat, nsims = 1000, 
                               lam = c(0.00117, 0.00137, 0.0045, 0.00193, 0.00174, 
                                       0.00106, 0.000783))
# save(sim_beans_prec, file = "saved_sims/sim_beans_prec.RData")

## Qda shrink sims Repeated 10 fold CV
sim_beans_qda <- qda_shrink_sim(dat, nsims = 1000, 
                                lam = c(0.00117, 0.00137, 0.0045, 0.00193, 0.00174, 
                                        0.00106, 0.000783))
# save(sim_beans_qda, file = "saved_sims/sim_beans_qda.RData")

## Median & se
load(file = "saved_sims/sim_beans_qda.RData")
sim_beans_qda |> 
  apply(2, function(x) paste0(round(median(x)*100, 2),"(", round(sd(x)*100, 3), ")"))

# HLDR Matrices -----------------------------------------------------------
eig_vals <- dat |> 
  group_by(class) |> 
  group_split() |> 
  lapply(function(x) eigen(solve(cov(x[,-1])))$value)

eig_vals |> 
  do.call(rbind, args = _) |> 
  as.data.frame() |> 
  pivot_longer(everything()) |> 
  mutate("eig_vect" = rep(1:16, times = 7), 
         "class" = rep(1:7, each = n()/7)) |> 
  ggplot(aes(eig_vect, value))+
  geom_point()+
  geom_line()+
  facet_wrap(~class, scales = "free_y")

lams <- c(eig_vals[[1]][2], eig_vals[[2]][4],  
          eig_vals[[3]][4], eig_vals[[4]][3], 
          eig_vals[[5]][2], eig_vals[[6]][2], 
          eig_vals[[7]][4])
## Saving HLDR Matrices for analysis
sim_beans_mats <- hldr_mats_sim(dat, 1000, lam =  lams, dims = c(15, 15, 7, 15, 5)) 
# save(sim_beans_mats, file = "saved_sims/sim_beans_mats.RData")

# HLDR Shrinkage ----------------------------------------------------------

## hldr shrink Repeated 10 fold CV
sim_beans_hldr <- hldr_sim(data = dat, 
                           lam = c(eig_vals[[1]][2], eig_vals[[2]][4],  
                                   eig_vals[[3]][4], eig_vals[[4]][3], 
                                   eig_vals[[5]][2], eig_vals[[6]][2], 
                                   eig_vals[[7]][4]), 
                           nsims = 1000)
# save(sim_beans_hldr, file = "saved_sims/sim_beans_hldr.RData")

load(file = "saved_sims/sim_beans_hldr.RData")
lapply(sim_beans_hldr, function(x) {
  meds <- apply(x, 2 , median)
  ses <- apply(x, 2, sd)
  min_dim <- which.min(meds)
  paste0(round(meds[min_dim], 4)*100, "(", round(ses[min_dim], 5)*100, ")", "(", min_dim, ")")
}) |> 
  unlist()

sim_hldr_boxplot(sim_beans_hldr, dims = 4:13)
load(file = "saved_sims/sim_beans_comp.RData")
sim_beans_comp <- comp_sim(data = dat, ends_u = 5, sir_d = 4, nsims = 1000)
apply(sim_beans_comp, 2, 
      function(y) paste0(round(median(y), 4), "(", round(sd(y), 4), ")"))
# save(sim_beans_comp, file = "saved_sims/sim_beans_comp.RData")
(apply(sim_beans_comp, 2, median)*100) |> round(2)

