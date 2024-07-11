library("tidyverse");theme_set(theme_minimal());library("qqplotr")
library("hldr")
library("SCPME")

# Getting and Cleaning Data -----------------------------------------------
set.seed(1)
dat <- iris |> 
  rename("class" = last_col()) |> 
  relocate(last_col(), .before = everything())
dat <- read_delim("sonar.csv", col_names = FALSE, delim = ",", na = "?") |> 
  janitor::clean_names() |> 
  rename("class" = last_col()) |> 
  relocate(class, .before = everything()) |>
  mutate(across(where(is.character), as_factor),
         across(-class, as.numeric)
         # "class" = case_when(class  == "0" ~ 1,
         #                     class =="1" ~ 2) |>
         #   as_factor()
  ) |>
  # select(-x14, -x15, -x28, -x16, -x21, -x22, -x27, -x26, -25) |> 
  drop_na()

out <- dat |> summary()
out[7,] |> as.list()
  
dat <- palmerpenguins::penguins |> 
  select(-year) |> 
  rename("class" = species) |> 
  mutate(across(where(is.character), as_factor), 
         across(-class, as.numeric)) |> 
  drop_na() 
  

dat$x17 <- dat$x17 + rnorm(length(dat$x17), sd = .00001)

dat[,c(18, 20)] <- dat[,c(18, 20)] + rnorm(length(dat[,c(18, 20)])*nrow(dat), sd = .001)

dat[,-1] <- dat[,-1] + rnorm(length(dat[,-1])^2, sd = .00001)

# Simulation --------------------------------------------------------------

dat |> 
  group_by(class) |> 
  count()
## Get lambda values for SCPME
x <- now()
omega <- SCPME_qda(dat[,-1], dat$class, nlam = 1000, cores = 7, K = 10)
lambdas <- vapply(seq_along(omega), function(i) omega[[i]]$Tuning[[2]], 
                  FUN.VALUE = numeric(1))
now() - x

# Repeated 10-fold cv
nsims <- 250
{
  x <- now()
  sim_div <- parallel::mclapply(1:nsims, function(sim){
    set.seed(sim)
    # 10-fold cv split
    splits <- (dat |> rsample::vfold_cv(strata = class))$splits
    train <- lapply(splits, rsample::training)
    test <- lapply(splits, rsample::testing)
    
    methods <- list(NULL, "Haff", "BGP16", "RidgeShrinkage", "SCPME")
    
    # SYS classifier model with training
    results <- lapply(seq_along(train), function(j){
      lapply(seq_along(methods), function(i){
        qda_shrink(train[[j]][,-1], train[[j]]$class, xtest = test[[j]][,-1], prec.est = methods[[i]], lam = lambdas)
      })
    })
    
   errors <- lapply(seq_along(methods), function(i){
      lapply(seq_along(results), function(k){
        mean(results[[k]][[i]] != test[[k]]$class)
      }) |> 
        do.call(rbind, args = _)
    }) |> 
      lapply(mean) |> 
      do.call(cbind, args = _)
   colnames(errors) <- c("QDA", methods[-1])
   errors
    
  }, mc.cores = parallel::detectCores() - 1) |> 
    do.call(rbind, args = _)
  print(now() - x)
}

sim_div |> 
  apply(2, median, na.rm = TRUE)

sim_div |> 
  as_tibble() |> 
  pivot_longer(everything()) |> 
  ggplot(aes(name, value, fill = name))+
  geom_boxplot()


# Dim Reduct --------------------------------------------------------------
nsims <- 250
{
  x <- now()
  sim_crx <- parallel::mclapply(1:nsims, function(sim){
    set.seed(sim)
    # 10-fold cv split
    splits <- (dat |> 
                 rsample::vfold_cv(strata = class))$splits
    train <- lapply(splits, rsample::training)
    test <- lapply(splits, rsample::testing)
    
    methods <- list(NULL, "Haff", "BGP16", "RidgeShrinkage", "SCPME")
    
    # SYS classifier model with training
    mods <- lapply(seq_along(train), function(j){
      lapply(seq_along(methods), function(i){
        hldr(class~., data = train[[j]], method = "SYS", prec.est = methods[[i]], lam = lambdas, fitted.class = FALSE)
      })
    })
    
    results <- lapply(seq_along(mods), function(j){
      lapply(seq_along(mods[[j]]), function(i){
        lapply(1:(ncol(train[[j]]) - 2), function(r){
          predict(mods[[j]][[i]], dims = 1:r,
                  newdata = test[[j]][,-1], grouping = test[[j]]$class)$class
        })
      })
    })
    
    for(j in seq_along(results)){
      results[[j]][[length(results[[j]]) + 1]] <- list(
        (MASS::qda(class~., data = train[[j]]) |>
           predict(newdata = test[[j]]))$class)
    }
    
    
    # cer by method and number of dims
    errors <- lapply(seq_along(mods), function(k){
      lapply(seq_along(results[[k]]), function(i){
        lapply(seq_along(results[[k]][[i]]), function(j){
          mean(results[[k]][[i]][[j]] != test[[k]]$class)
        }) |> 
          do.call(cbind, args = _)
      }) 
    })
    
    # rbind cer by dim for each fold and compute mean by each dim
    lapply(seq_along(errors[[1]]), function(i) {
      lapply(seq_along(errors), function(j) {
        errors[[j]][[i]]
      }) |>
        do.call(rbind, args = _)
    }) |>
      lapply(colMeans)
  }, mc.cores = parallel::detectCores() - 1)
  print(now() - x)
  # rbind cer by method
  sim_crx <- lapply(seq_along(sim_crx[[1]]), function(i){
      lapply(seq_along(sim_crx), function(j){
        sim_crx[[j]][[i]]
      }) |> 
        do.call(rbind, args = _)
    })
  names(sim_crx) <- c("none", "Haff", "BGP16", "RidgeShrinkage", "SCPME", "QDA")
}
lapply(sim_crx, function(x) {
  apply(x, 2, function(y) sd(y)) |>
    round(5) |> 
    min()
})











####### Sonar







# Boxplots ----------------------------------------------------------------
nms <-  names(sim_crx[which(unlist(lapply(sim_crx, ncol)) != 1)])
sim_crx_long <- sim_crx[which(unlist(lapply(sim_crx, ncol)) != 1)] |> 
  do.call(rbind, args = _) |> 
  as_tibble(.name_repair  = "minimal") |> 
  set_names(paste0(1:ncol(sim_crx[[1]]))) |> 
  mutate("method" = rep(nms, each = n()/length(nms)) |> 
           as_factor()) |> 
  pivot_longer(cols = -method,
               names_to = "ndims",
               values_to = "cer") |>
  mutate("ndims" = as_factor(ndims)) |> 
  filter(as.numeric(ndims) %in% 8:30)
# df of outliers
outliers <- sim_crx_long |>
  group_by(method, ndims) |>
  filter(cer > 1.5*IQR(cer) + quantile(cer, .75) |
           cer < 1.5*IQR(cer) - quantile(cer, .75)) |>
  ungroup() |>
  mutate("x_dum" = case_when(method == "none" ~ as.numeric(ndims) - .28,
                             method == "Haff" ~ as.numeric(ndims) - .1,
                             method == "BGP16" ~ as.numeric(ndims) + .1,
                             method == "RidgeShrinkage" ~ as.numeric(ndims) + .28, 
                             method == "SCPME" ~ as.numeric(ndims) + .38
                             )
  )


meds <- sim_crx[which(unlist(lapply(sim_crx, ncol)) == 1)] |> 
  lapply(function(x) median(x)) |> 
  do.call(bind_cols, args = _) |> 
  pivot_longer(everything(), names_to = "method") |> 
  arrange(method) |> 
  mutate("method" = as_factor(method)) |> 
  filter(method != "LDA")
# boxplot by method and ndims
sim_crx_long |> 
  ggplot(aes(ndims, cer, fill = method, shape = method))+
  stat_boxplot(geom='errorbar')+
  geom_boxplot(outlier.shape = NA)+
  # geom_point(aes(x = x_dum), data = outliers)+
  geom_hline(aes(yintercept = value, color = method, linetype = method), 
             data = meds)+
  scale_fill_manual(values = c("purple", "#31688e", "#35b779", "#fde725", "blue"))+
  scale_shape_manual(values = c(16, 17, 15, 18, 19))+
  scale_color_manual(values = c("red"), 
                     labels = meds$method)+
  scale_linetype_manual(values = c(1), 
                        labels =  meds$method)+
  labs(x = "Number of Dimensions",
       y = "CER",
       color = "", linetype = "",
       fill = "Method", shape = "Method",
       title = "Figure ?: Simulation ? CER Boxplots by Dimension Ranking")

