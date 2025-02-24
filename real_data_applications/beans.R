source("wrapper_fns.R")

#  Reading Data & Cleaning ------------------------------------------------
dat <- read_csv("datasets/Dry_Bean_Dataset.arff", col_names = FALSE) |> 
  janitor::clean_names() |>
  rename("class" = last_col()) |> 
  relocate(class, .before = everything()) |>
  mutate(across(where(is.numeric), function(x)(x - mean(x))/sd(x)), 
         "class" = as_factor(class) |> as.numeric() |> as_factor())

# Comparing Prec Est in SDRS ----------------------------------------------------------
## Generating tuning params
eig_vals <- dat |> 
  group_by(class) |> 
  group_split() |> 
  lapply(function(x) eigen(solve(cov(x[,-1])))$value)

lams <- c(eig_vals[[1]][2], eig_vals[[2]][4],  
          eig_vals[[3]][4], eig_vals[[4]][3], 
          eig_vals[[5]][2], eig_vals[[6]][2], 
          eig_vals[[7]][4])

## SDRS comparison
sim_beans_sdrs <- sdrs_sim(data = dat, lambdas = lams, nsims = 1000)
# save(sim_beans_sdrs, file = "saved_sims/sim_beans_sdrs.RData")

## Median & SE
load(file = "saved_sims/sim_beans_sdrs.RData")
lapply(sim_beans_sdrs, function(x) {
  meds <- apply(x, 2 , median)
  ses <- apply(x, 2, sd)
  min_dim <- which.min(meds)
  paste0(round(meds[min_dim], 4)*100, "(", round(ses[min_dim], 5)*100, ")", "(", min_dim, ")")
}) |> 
  unlist()

## Boxplots
sim_sdrs_boxplot(sim_beans_sdrs, dims = 4:13)


# Competitors -------------------------------------------------------------
sim_beans_comp <- comp_sim(data = dat, ends_u = 5, sir_d = 4, nsims = 1000)
# save(sim_beans_comp, file = "saved_sims/sim_beans_comp.RData")
load(file = "saved_sims/sim_beans_comp.RData")
apply(sim_beans_comp, 2, 
      function(y) paste0(round(median(y), 4), "(", round(sd(y), 4), ")"))


