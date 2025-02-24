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

## adding error for singularity issues
set.seed(1)
dat[,-1] <- dat[,-1] + matrix(rnorm(nrow(dat[,-1])*ncol(dat[,-1]), sd = .00001), nrow = nrow(dat))


# Comparing Prec Est in SDRS ----------------------------------------------------------
sim_div_sdrs <- sdrs_sim(data = dat, lam = c(1.5, 1.5), nsims = 1000, 
                         type = "qda", gamma = c(.55, 3.04))
# save(sim_div_sdrs, file = "saved_sims/sim_div_sdrs.RData")

## Median & SE
load(file = "saved_sims/sim_div_sdrs.RData")
lapply(sim_div_sdrs, function(x) {
  meds <- apply(x, 2 , median)
  ses <- apply(x, 2, sd)
  min_dim <- which.min(meds)
  paste0(round(meds[min_dim], 4)*100, "(", round(ses[min_dim], 5)*100, ")", "(", min_dim, ")")
}) |> 
  unlist()

## Boxplots
sim_sdrs_boxplot(sim_div_sdrs, dims = 2:9)

# Competitors -------------------------------------------------------------
sim_div_comp <- comp_sim(dat, ends_u = 1, nsims = 1000)
# save(sim_div_comp, file = "saved_sims/sim_div_comp.RData")
load(file = "saved_sims/sim_div_comp.RData")
apply(sim_div_comp, 2, 
      function(y) paste0(round(median(y), 4), "(", round(sd(y), 4), ")"))
