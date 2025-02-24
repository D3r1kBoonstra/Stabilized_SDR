source("wrapper_fns.R")

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

# Comparing Prec Est in SDRS ----------------------------------------------------------
sim_ion_sdrs <- sdrs_sim(data = dat, lam = c(.0005, .0005), nsims = 1000, 
                        gamma = c(1.04, .43), type = "qda", 
                        standardize_xbar = FALSE)
# save(sim_ion_sdrs, file = "saved_sims/sim_ion_sdrs.RData")

## Median & SE
load(file = "saved_sims/sim_ion_sdrs.RData")
lapply(sim_ion_sdrs, function(x) {
  meds <- apply(x, 2 , median)
  ses <- apply(x, 2, sd)
  min_dim <- which.min(meds)
  paste0(round(meds[min_dim], 5)*100, "(", round(ses[min_dim], 5)*100, ")", "(", min_dim, ")")
}) |> 
  unlist()

## Boxplots
sim_sdrs_boxplot(sim_ion_sdrs, dims = 2:11)


# Competitors -------------------------------------------------------------
sim_ion_comp <- comp_sim(data = dat, ends_u = 1, nsims = 1000)
# save(sim_ion_comp, file = "saved_sims/sim_ion_comp.RData")
load(file = "saved_sims/sim_ion_comp.RData")
apply(sim_ion_comp, 2, 
      function(y) paste0(round(median(y), 4), "(", round(sd(y), 4), ")"))

