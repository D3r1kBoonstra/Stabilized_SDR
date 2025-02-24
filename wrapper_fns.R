library("tidyverse");theme_set(theme_minimal())
# MC Sims -----------------------------------------------------------------
## Wrapper fn to perform MC sim in paper
sim_wrapper <- function(mu, Sigma, n, ni, nsims, lambdas, type = "L1", gamma = NULL){
  p <- length(mu[[1]])
  
  x <- now()
  out <- parallel::mclapply(1:nsims, function(sim){
    set.seed(sim)
    # generating sample
    data <- lapply(seq_along(mu), function(i) {
      cbind("class" = i, MASS::mvrnorm(n, mu[[i]], Sigma[[i]])) |>
        as.data.frame()
    }) |> 
      do.call(rbind, args = _)
    
    lapply(ni, function(ni){
      set.seed(sim)
      splits <- data |> 
        rsample::initial_split(strata = class, prop = ni/n)
      train <- rsample::training(splits)
      test <- rsample::testing(splits)
      
      methods <- list(NULL, "Haff", "Wang", "Bodnar", "MRY")
      
      mods <- lapply(seq_along(methods), function(i){
        sdr::sdr(as.matrix(train[,-1]), grouping = train$class, method = "SDRS", prec.est = methods[[i]], 
             lam = lambdas, type = type, gamma = gamma)
      })
      
      
      qdf_err <- mean((MASS::qda(class~., data = train) |>
                         predict(newdata = test))$class != test$class)
      
      lapply(seq_along(mods), function(i){
        sapply(1:p, function(r){
          pred <- predict(mods[[i]], dims = 1:r,
                          newdata = test[,-1], grouping = test$class)$class
          mean(pred != test$class)
        }) 
      }) |> 
        do.call(rbind, args = _) |> 
        as.data.frame() |> 
        rbind(rep(qdf_err, times = p)) |> 
        mutate(c("MLE", "Haff", "Wang", "Bodnar", "MRY", "qdf"), 
               rep(ni, times = length(methods) + 1)
        )
    }) |> 
      do.call(rbind, args = _)
  }, mc.cores = parallel::detectCores() - 1) |> 
    do.call(rbind, args = _) |> 
    setNames(c(paste0(1:(p)), "method", "ni"))
  message(now() - x)
  out
}

## Summarizes MC Sim by CER, SE and min dim
sim_summary <- function(sim){
  sim |> 
    mutate("method" = as_factor(method)) |> 
    group_by(ni) |> 
    group_split() |> 
    lapply(function(x){
      x <- x |> 
        group_by(method) |> 
        group_split() |> 
        setNames(c("S_inv", "Haff", "Wang", "Bodnar", "MRY", "QDA")) |> 
        lapply(function(x) {
          x <- x |> 
            select(-method, -ni)
          meds <- apply(x, 2 , median)
          ses <- apply(x, 2, sd)
          min_dim <- which.min(meds)
          paste0(round(meds[min_dim], 4)*100, "(", round(ses[min_dim], 5)*100, ")", "(", min_dim, ")")
        }) |> 
        unlist()
    }) |> 
    setNames(paste0("ni = ", unique(sim$ni)))
}

## Plots MC Sim as boxplots by dim
sim_boxplots <- function(sim){
  qdf_meds <- sim |> 
    filter(method == "qdf") |> 
    summarise("med" = median(`1`), .by = ni)
  
  sim |> 
    filter(method != "qdf") |> 
    pivot_longer(cols = -c(method, ni), 
                 names_to = "dim", 
                 values_to = "error") |> 
    mutate("method" = as_factor(method), 
           "dim" = as_factor(dim)) |> 
    ggplot(aes(dim, error, fill =  method))+
    geom_boxplot(outliers = FALSE)+
    geom_hline(aes(yintercept = med, color = "1"), data = qdf_meds)+
    scale_color_manual(values = "red", labels = "QDA")+
    facet_wrap(~ni, nrow = 3, scales = "free_y")+
    scale_fill_manual(values = c("purple", "#6495ED", "#21918c", "#5ec962", "#fde725"))+
    labs(y = expression(bar("CER")), 
         color = "", fill = "Precision", shape = "Precision")+
    theme(strip.background =element_rect(fill="lightgrey"), 
          legend.position = "bottom")
}

# Real Data Application ----------------------------------------------------------
## Wrapper fn to perform real data applications for comparing prec est in SDRS
sdrs_sim <- function(data, lambdas, nsims, gamma = NULL, type = "L1", standardize_xbar = FALSE){
  x <- now()
  sim <- parallel::mclapply(1:nsims, function(sim){
    set.seed(sim)
    # 10-fold cv split
    splits <- (data |> rsample::vfold_cv(strata = class))$splits
    train <- lapply(splits, rsample::training)
    test <- lapply(splits, rsample::testing)
    
    methods <- list(NULL, "Haff", "Wang", "Bodnar" , "MRY")
    
    # SYS classifier model with training
    mods <- lapply(seq_along(train), function(j){
      lapply(seq_along(methods), function(i){
        sdr::sdr(train[[j]][,-1], grouping = train[[j]]$class, method = "SDRS", prec.est = methods[[i]], 
            lam = lambdas, gamma = gamma, 
            type = type, standardize_xbar = standardize_xbar)
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
  message(paste0("Time Difference of ", now() - x))
  # rbind cer by method
  sim <- lapply(seq_along(sim[[1]]), function(i){
    lapply(seq_along(sim), function(j){
      sim[[j]][[i]]
    }) |> 
      do.call(rbind, args = _)
  })
  names(sim) <- c("none", "Haff", "Wang", "Bodnar", "Molstad", "QDA")
  sim
}

## Boxplot for comparing prec ests in SDRS
sim_sdrs_boxplot <- function(sim, dims = NULL){
  if(is.null(dims)) dims <- 1:(lapply(sim, ncol) |> unlist() |> max())
  nms <-  names(sim[which(unlist(lapply(sim, ncol)) != 1)])
  sim_long <- sim[which(unlist(lapply(sim, ncol)) != 1)] |> 
    do.call(rbind, args = _) |> 
    as_tibble(.name_repair  = "minimal") |> 
    set_names(paste0(1:ncol(sim[[1]]))) |> 
    mutate("method" = rep(nms, each = n()/length(nms)) |> 
             as_factor() |> 
             fct_recode("MLE" = "none", "MRY" = "Molstad")) |> 
    pivot_longer(cols = -method,
                 names_to = "ndims",
                 values_to = "cer") |>
    filter(as.numeric(ndims) %in% dims) |> 
    mutate("ndims" = as_factor(ndims)) 
  
  # df of outliers
  outliers <- sim_long |>
    mutate("x_dum" = case_when(method == "MLE" ~ as.numeric(ndims) - .3,
                               method == "Haff" ~ as.numeric(ndims) - .15,
                               method == "Wang" ~ as.numeric(ndims),
                               method == "Bodnar" ~ as.numeric(ndims) + .15,
                               method == "MRY" ~ as.numeric(ndims) + .3)
    )|> 
    group_by(method, ndims) |>
    filter(cer > 1.5*IQR(cer) + quantile(cer, .75) |
             cer < 1.5*IQR(cer) - quantile(cer, .75)) |>
    ungroup()
  
  
  meds <- sim[which(unlist(lapply(sim, ncol)) == 1)] |> 
    lapply(function(x) median(x)) |> 
    do.call(bind_cols, args = _) |>
    rename("QDF" = "QDA") |> 
    pivot_longer(everything(), names_to = "method") |> 
    mutate("method" = as_factor(method))
  # boxplot by method and ndims
  sim_long |>
    ggplot(aes(ndims, cer, fill = method, shape = method))+
    stat_boxplot(geom='errorbar')+
    geom_boxplot(outlier.shape = NA, linewidth = .2)+
    geom_point(aes(x = x_dum, cer, shape = method), data = outliers, 
               alpha = .5)+
    geom_hline(aes(yintercept = value, color = "1"), data = meds)+
    scale_color_manual(values = "red", labels = "QDF")+
    scale_fill_manual(values = c("purple", "#6495ED", "#21918c", "#5ec962", "#fde725"))+
    scale_shape_manual(values = c(16, 17, 15, 3, 7))+
    labs(x = "Dimension (r)",
         y = expression(bar("CER")),
         color = "", linetype = "",
         fill = "Precision", shape = "Precision")+
    theme(legend.position = "bottom")+
    guides(color = guide_legend(order = 1), 
           shape = guide_legend(override.aes = list(alpha = 1)))
  
}

## Updating Lasso Sir fn to fix bugs
lassosir_fixed <- function (X, Y, H = 0, choosing.d = "automatic", solution.path = FALSE, 
                            categorical = FALSE, nfolds = 10, screening = TRUE, no.dim = 0) 
{
  if (no.dim != 0) 
    choosing.d = "given"
  if ((categorical == FALSE) & (H == 0)) {
    H <- (function() {
      H <- readline("For the continuous response, please choose the number of slices:   ")
      H <- as.numeric(unlist(strsplit(H, ",")))
      return(dim)
    })()
  }
  p <- dim(X)[2]
  n <- dim(X)[1]
  if (categorical == FALSE) {
    ORD <- order(Y)
    X <- X[ORD, ]
    Y <- Y[ORD]
    ms <- array(0, n)
    m <- floor(n/H)
    c <- n%%H
    M <- matrix(0, nrow = H, ncol = n)
    if (c == 0) {
      M <- diag(H) %x% matrix(1, nrow = 1, ncol = m)/m
      ms <- m + ms
    }
    else {
      for (i in 1:c) {
        M[i, ((m + 1) * (i - 1) + 1):((m + 1) * i)] <- 1/(m + 
                                                            1)
        ms[((m + 1) * (i - 1) + 1):((m + 1) * i)] <- m
      }
      for (i in (c + 1):H) {
        M[i, ((m + 1) * c + (i - c - 1) * m + 1):((m + 
                                                     1) * c + (i - c) * m)] <- 1/m
        ms[((m + 1) * c + (i - c - 1) * m + 1):((m + 
                                                   1) * c + (i - c) * m)] <- m - 1
      }
    }
    if (screening == TRUE) {
      x.sliced.mean <- M %*% X
      sliced.variance <- apply(x.sliced.mean, 2, var)
      keep.ind <- sort(order(sliced.variance, decreasing = TRUE)[1:n])
    }
    else {
      keep.ind <- c(1:p)
    }
    X <- X[, keep.ind]
    X.H <- matrix(0, nrow = H, ncol = dim(X)[2])
    grand.mean <- matrix(apply(X, 2, mean), nrow = 1, ncol = dim(X)[2])
    X.stand.ord <- X - grand.mean %x% matrix(1, nrow = dim(X)[1], 
                                             ncol = 1)
    X.H <- M %*% X.stand.ord
  }
  else {
    ms <- array(0, n)
    Y.unique <- unique(Y)
    H <- length(Y.unique)
    ORD <- which(Y == Y.unique[1])
    nH <- sum(Y == Y.unique[1])
    ms[1:nH] <- nH
    for (i in 2:H) {
      ORD <- c(ORD, which(Y == Y.unique[i]))
      nH <- c(nH, sum(Y == Y.unique[i]))
      ms[(sum(nH[1:(i - 1)]) + 1):sum(nH[1:i])] <- nH[i]
    }
    X <- X[ORD, ]
    M <- matrix(0, nrow = H, ncol = n)
    M[1, 1:nH[1]] <- 1/nH[1]
    for (i in 2:H) M[i, (sum(nH[1:(i - 1)]) + 1):sum(nH[1:i])] <- 1/nH[i]
    if (screening == TRUE) {
      x.sliced.mean <- M %*% X
      sliced.variance <- apply(x.sliced.mean, 2, var)
      keep.ind <- sort(order(sliced.variance, decreasing = TRUE)[1:n])
    }
    else {
      keep.ind <- c(1:p)
    }
    X <- X[, keep.ind]
    X.H <- matrix(0, nrow = H, ncol = dim(X)[2])
    grand.mean <- matrix(apply(X, 2, mean), nrow = 1, ncol = dim(X)[2])
    X.stand.ord <- X - grand.mean %x% matrix(1, nrow = dim(X)[1], 
                                             ncol = 1)
    X.H <- M %*% X.stand.ord
  }
  svd.XH <- svd(X.H, nv = p)
  res.eigen.value <- array(0, p)
  res.eigen.value[1:dim(X.H)[1]] <- (svd.XH$d)^2/H
  if (choosing.d == "manual") {
    plot(c(1:p), res.eigen.value, ylab = "eigen values")
    no.dim <- (function() {
      dim <- readline("Choose the number of directions:   ")
      dim <- as.numeric(unlist(strsplit(dim, ",")))
      return(dim)
    })()
  }
  if (choosing.d == "automatic") {
    beta.hat <- array(0, c(p, min(p, H)))
    Y.tilde <- array(0, c(n, min(p, H)))
    for (ii in 1:min(p, H)) {
      eii <- matrix(0, nrow = dim(svd.XH$v)[2], ncol = 1)
      eii[ii] <- 1
      eigen.vec <- solve(t(svd.XH$v), eii)
      Y.tilde[, ii] <- t(M) %*% M %*% X.stand.ord %*% 
        eigen.vec/(res.eigen.value[ii]) * matrix(1/ms, 
                                                 nrow = n, ncol = 1)
    }
    mus <- array(0, min(p, H))
    for (ii in 1:min(p, H)) {
      lars.fit.cv <- cv.glmnet(X.stand.ord, Y.tilde[, 
                                                    ii], nfolds = nfolds)
      ind <- max(which(lars.fit.cv$cvm == min(lars.fit.cv$cvm)))
      if (ind == 1) 
        ind <- 2
      lambda <- lars.fit.cv$lambda[ind]
      mus[ii] <- lambda
      lars.fit <- glmnet(X.stand.ord, Y.tilde[, ii], lambda = lambda)
      beta.hat[keep.ind, ii] <- as.double(lars.fit$beta)
    }
    temp.2 <- sqrt(apply(beta.hat^2, 2, sum)) * res.eigen.value[1:H]
    temp <- temp.2/temp.2[1]
    res.kmeans <- kmeans(temp, centers = 2, algorithm = "Lloyd")
    no.dim <- min(sum(res.kmeans$cluster == 1), sum(res.kmeans$cluster == 
                                                      2))
  }
  beta.hat <- array(0, c(p, no.dim))
  Y.tilde <- array(0, c(n, no.dim))
  for (ii in 1:no.dim) {
    eii <- matrix(0, nrow = dim(t(svd.XH$v))[2], ncol = 1)
    eii[ii] <- 1
    eigen.vec <- solve(t(svd.XH$v), eii)
    Y.tilde[, ii] <- t(M) %*% M %*% X.stand.ord %*% eigen.vec/(res.eigen.value[ii]) * 
      matrix(1/ms, nrow = n, ncol = 1)
  }
  if (solution.path == FALSE) {
    mus <- array(0, no.dim)
    if (no.dim == 1) {
      lars.fit.cv <- cv.glmnet(X.stand.ord, Y.tilde, nfolds = nfolds)
      ind <- max(which(lars.fit.cv$cvm == min(lars.fit.cv$cvm)))
      if (ind == 1) 
        ind <- 2
      lambda <- lars.fit.cv$lambda[ind]
      lars.fit <- glmnet(X.stand.ord, Y.tilde, lambda = lambda)
      beta.hat[keep.ind] <- as.double(lars.fit$beta)
    }
    else {
      for (ii in 1:no.dim) {
        lars.fit.cv <- cv.glmnet(X.stand.ord, Y.tilde[, 
                                                      ii], nfolds = nfolds)
        ind <- max(which(lars.fit.cv$cvm == min(lars.fit.cv$cvm)))
        if (ind == 1) 
          ind <- 2
        lambda <- lars.fit.cv$lambda[ind]
        mus[ii] <- lambda
        lars.fit <- glmnet(X.stand.ord, Y.tilde[, ii], 
                           lambda = lambda)
        beta.hat[keep.ind, ii] <- as.double(lars.fit$beta)
      }
    }
    list(beta = beta.hat, eigen.value = res.eigen.value, 
         no.dim = no.dim, H = H, categorical = categorical)
  }
  else {
    lars.fit.all <- list()
    for (ii in 1:no.dim) {
      lars.fit.all[[ii]] <- glmnet(X.stand.ord, Y.tilde[, 
                                                        ii])
    }
    lars.fit.all
  }
}
environment(lassosir_fixed) <- asNamespace('LassoSIR')

## Wrapper fn to compared SDRS to competive methods for Real data
comp_sim <- function(data, ends_u = 1, sir_d = 0, nsims, ncores = parallel::detectCores() - 1){
  x <- now()
  filedr <- rstudioapi::getActiveDocumentContext()$path |>
    str_remove("/[^/]+$")
  sim_data_pth <- file.path(filedr, "ends_sim")
  dir.create(sim_data_pth, showWarnings = FALSE)
  unlink(file.path(sim_data_pth, "*"), recursive = TRUE, force = TRUE)
  write_csv(as.data.frame(ends_u), file.path(sim_data_pth, "u.csv"))
  write_csv(as.data.frame(ncores), file.path(sim_data_pth, "ncores.csv"))
  write_csv(as.data.frame(nsims), file.path(sim_data_pth, "nsims.csv"))
  write_csv(data, file.path(sim_data_pth, "data.csv"))
  
 
  
  errors <- parallel::mclapply(1:nsims, function(j){
    set.seed(j)
    splits <- (data |> rsample::vfold_cv(strata = class))$splits
    train <- lapply(splits, rsample::training)
    test <- lapply(splits, rsample::testing)
    
    lapply(seq_along(train), function(i){
      
      sir.lasso <- lassosir_fixed(as.matrix(train[[i]][,-1]), train[[i]]$class, choosing.d="automatic",
                                  solution.path=FALSE, categorical = TRUE, nfolds=10,
                                  screening=FALSE, no.dim = sir_d)
      sir_pred <- (MASS::qda(x = sdr::project(as.matrix(train[[i]][,-1]), t(sir.lasso$beta)), grouping = train[[i]]$class) |> 
                     predict(newdata = sdr::project(as.matrix(test[[i]][,-1]), t(sir.lasso$beta))))$class
      
      msda_cv <- TULIP::cv.msda(as.matrix(train[[i]][,-1]), y = as.numeric(train[[i]]$class), nfolds = 10)
      msda_pred <- TULIP::msda(as.matrix(train[[i]][,-1]), y = as.numeric(train[[i]]$class), testx = as.matrix(test[[i]][,-1]), lambda = msda_cv$lambda.min)$pred
      
      if(length(unique(data$class)) == 2){
        qdap_pred <- QDAP::qdap(as.matrix(train[[i]][,-1]), y = ifelse(train[[i]]$class == 1, 1, 0), xnew = test[[i]][,-1])$class
      } 
      
      c(
        "msda" = mean(as.numeric(test[[i]]$class) != msda_pred), 
        "qdap" = if(length(unique(data$class)) == 2){
          mean(ifelse(test[[i]]$class == 1, 1, 0) != qdap_pred)
        } else {
          NA
        }, 
        "lassoSir" = mean(test[[i]]$class != sir_pred)
      )
    }) |> 
      do.call(rbind, args = _) |> 
      colMeans()
    
  }, mc.cores = ncores) |> 
    do.call(rbind, args = _)
  
  message("starting ENDS Matlab")
  
  script_path <- here::here("real_data_applications", "ends.m")
  
  ## Change Matlab location to yours
  matlab_cmd <- sprintf(
    '/Applications/MATLAB_R2023b.app/bin/matlab -nodisplay -batch "run(\'%s\'); exit"',
    script_path
  )
  # Run MATLAB from R
  system(matlab_cmd)
  
  ends_error <- read_csv(file.path(sim_data_pth, "ENDS_error.csv"),
                         col_names = FALSE,
                         show_col_types = FALSE)[[1]]
  
  message(now() - x)
  cbind(errors, "ends" = ends_error)
}