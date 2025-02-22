library("tidyverse");library("qqplotr");theme_set(theme_minimal())
library("sdr")
# MC Sims -----------------------------------------------------------------
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
      
      methods <- list(NULL, "Haff", "RidgeShrinkage", "BGP16", "SCPME")
      
      mods <- lapply(seq_along(methods), function(i){
        sdr(as.matrix(train[,-1]), grouping = train$class, method = "SDRS", prec.est = methods[[i]], 
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

prec_est_sim2 <- function(mu, Sigma, n, ni, nsims, lambdas, gamma = NULL, type = "L1", standardize_xbar = FALSE){
  x <- now()
  sim <- parallel::mclapply(1:nsims, function(sim){
    set.seed(sim)
    
    data <- lapply(seq_along(mu), function(i) {
      cbind("class" = i, MASS::mvrnorm(n, mu[[i]], Sigma[[i]])) |>
        as.data.frame()
    }) |> 
      do.call(rbind, args = _)
    
    lapply(ni, function(ni){
      
      splits <- data |> 
        rsample::initial_split(strata = class, prop = ni/n)
      train <- rsample::training(splits)
      test <- rsample::testing(splits)
      
      methods <- list(NULL, "Haff", "RidgeShrinkage", "BGP16" , "SCPME")
      
      out <- lapply(seq_along(methods), function(i){
        prec_est(train[,-1], train$class, prec.est = methods[[i]],
                 lam = lambdas, gamma = gamma,
                 type = type, standardize_xbar = standardize_xbar)
      })
      # 
      # out <- lapply(seq_along(methods), function(i){
      #   sdr(class~., data = train, prec.est = methods[[i]], 
      #            lam = lambdas, gamma = gamma, 
      #            type = type, standardize_xbar = standardize_xbar)$ProjectionMatrix |> 
      #     as.matrix()
      # })
      
      # out[[5]]$ProjectionMatrix
      # 
      # dim(out[[5]]$M)
      # 
      # S_inv <- lapply(Sigma, solve)
      # M <- cbind(S_inv[[2]] %*% as.matrix(mu[[2]]) - S_inv[[1]] %*% as.matrix(mu[[1]]), 
      #            # S_inv[[3]] %*% as.matrix(mu[[3]]) - S_inv[[1]] %*% as.matrix(mu[[1]]),
      #            Sigma[[2]] - Sigma[[1]]
      #            # , Sigma[[3]] - Sigma[[1]]
      # )
      
      # out[[length(out) + 1]] <-  t(svd(M)$u) |> as.matrix()
      out[[length(out) + 1]] <- lapply(Sigma, solve)
      out
    })
    
  }, mc.cores = parallel::detectCores() - 1) 
  
  out <- lapply(seq_along(sim[[1]]), function(i){
    lapply(seq_along(sim[[1]][[1]]), function(j){
      lapply(seq_along(sim[[1]][[1]][[1]]), function(k){
        (lapply(seq_along(sim), function(l){
          sim[[l]][[i]][[j]][[k]]
        }) |> 
          Reduce("+", x = _))/length(sim)
      })
    })
  })
  # names(out) <- paste0("ni = ", ni)
  # for(i in seq_along(out)){
  #   names(out[[i]]) <- c("MLE", "Haff", "Wang", "Bodnar", "MRY", "TRUE")
  # }
  message(paste0("Time Difference of ", now() - x))
  system("say Just finished!")
  out
}

sim_summary <- function(sim){
  sim |> 
    mutate("method" = as_factor(method)) |> 
    group_by(ni) |> 
    group_split() |> 
    lapply(function(x){
      x <- x |> 
        group_by(method) |> 
        group_split() |> 
        setNames(c("MLE", "Haff", "Wang", "Bodnar", "MRY", "qdf")) |> 
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
    scale_color_manual(values = "red", labels = "QDF")+
    facet_wrap(~ni, nrow = 3, scales = "free_y")+
    scale_fill_manual(values = c("purple", "#6495ED", "#21918c", "#5ec962", "#fde725"))+
    labs(y = expression(bar("CER")), 
         color = "", fill = "Precision", shape = "Precision")+
    theme(strip.background =element_rect(fill="lightgrey"), 
          legend.position = "bottom")
}

# Precision Estimation ----------------------------------------------------
prec_est <- function(x, grouping, prec.est = NULL, ...) {
  x <- as.matrix(x)
  classes <- unique(grouping)
  data <- lapply(1:length(classes), function(i) {
    x[grouping == classes[[i]],]
  })
  p <- ncol(x)
  n <- lapply(data, nrow)
  S <- lapply(data, stats::cov)
  if(is.null(prec.est)) {
    S_inv <- lapply(S, function(x) solve(x, tol = NULL))
  } else if(prec.est == "Haff"){
    S_inv <- lapply(seq_along(S), function(i, ...){
      Haffshrinkage(S = S[[i]], n = n[[i]], tol = NULL)
    })
  } else if(prec.est == "RidgeShrinkage"){
    S_inv <- lapply(seq_along(S), function(i, ...){
      RidgeShrinkage(S = S[[i]], n = n[[i]], tol = NULL)
    })
  } else if(prec.est == "BGP16"){
    S_inv <- lapply(seq_along(S), function(i, ...){
      p <- ncol(S[[i]])
      HDShOP::InvCovShrinkBGP16(n[[i]], p, diag(p), solve(S[[i]], tol = NULL))$S
    })
  } else if(prec.est == "SCPME"){
    S_inv <- SCPME_qda(x, grouping, ...)
    S_inv <- lapply(seq_along(S_inv), function(i) S_inv[[i]]$Omega)
  }
  ## Return
  S_inv
}



prec_est_sim <- function(data, lambdas, nsims, gamma = NULL, type = "L1", standardize_xbar = FALSE){
  x <- now()
  sim <- parallel::mclapply(1:nsims, function(sim){
    set.seed(sim)
    # 10-fold cv split
    splits <- (data |> rsample::vfold_cv(strata = class))$splits
    train <- lapply(splits, rsample::training)
    test <- lapply(splits, rsample::testing)
    
    methods <- list(NULL, "Haff", "RidgeShrinkage", "BGP16" , "SCPME")
    
    # SYS classifier model with training
    results <- lapply(seq_along(train), function(j){
      lapply(seq_along(methods), function(i){
        prec_est(train[[j]][,-1], train[[j]]$class, prec.est = methods[[i]], 
                 lam = lambdas, gamma = gamma, 
                 type = type, standardize_xbar = standardize_xbar)
      })
    })
    
    lapply(seq_along(methods), function(i){
        lapply(seq_along(results[[1]][[1]]), function(k){
        (lapply(seq_along(results), function(j){
          results[[j]][[i]][[k]]
        }) |> 
          Reduce("+", x = _))/10
      })
    })
  }, mc.cores = parallel::detectCores() - 1) 
  out <- lapply(1:5, function(i){
    lapply(seq_along(sim[[1]][[1]]), function(k){
      ((lapply(seq_along(sim), function(j){
        sim[[j]][[i]][[k]]
      }) |> 
        Reduce("+", x = _))/nsims) |> 
        zapsmall()
    })
  })
  names(out) <- c("MLE", "Haff", "Wang", "Bodnar", "MRY")
  message(paste0("Time Difference of ", now() - x))
  system("say Just finished!")
  out
}

prec_analysis <- function(prec_list, title = "Dataset"){
  
  ## Helper fns
  upper_tri <- function(x){
    x[lower.tri(x)]<- NA
    #return
    x
  }
  ggarrange_nrow1 <- function(...){
    ggpubr::ggarrange(..., nrow = 1, align = "hv", widths = c(1, 1))
  }
  ggarrange_ncol1 <- function(...){
    ggpubr::ggarrange(..., ncol = 1, align = "v", heights = c(1, 1.25))
  }
  
  methods <- c("MLE", "Haff", "Wang", "Bodnar", "MRY", "TRUE")
  n_methods <- length(methods)
  k <- length(prec_list[[1]])
  
  eigs <- lapply(1:n_methods, function(i) {
    lapply(1:k, function(j) {
      eigen(prec_list[[i]][[j]])
    })
  })
  
  eig_vals <- lapply(1:n_methods, function(i){
    lapply(1:k, function(j){
      eigs[[i]][[j]]$values
    })|> 
      do.call(cbind, args = _) |> 
      as.data.frame() |> 
      setNames(paste0("Class", 1:k)) |> 
      mutate("vector" = 1:n())
  }) |> 
    do.call(rbind, args = _) |> 
    mutate("precision" = rep(methods,
                             each = n()/n_methods) |> 
             as_factor())
  
  eig_vects <- lapply(1:n_methods, function(i){
    lapply(1:k, function(j){
      eigs[[i]][[j]]$vectors
    }) 
  }) 
  
  plots <- lapply(seq_along(prec_list), function(i){
    temp <- lapply(seq_along(prec_list[[i]]), function(j){
      colors <- cbind(palette.colors(n = j, palette = "Set 1"),
                      palette.colors(n = j, palette = "Set 2"))
      # colors <- cbind(c("#FC766AFF", "#5F4B8BFF"),
      #                 c("#5B84B1FF", "forestgreen"))
      dimnames(prec_list[[i]][[j]]) <- NULL
      df <- prec_list[[i]][[j]] |> 
        round(5) |> 
        zapsmall() |> 
        upper_tri() |>
        reshape2::melt(na.rm = TRUE) |>
        mutate("class" = rep(paste0("Class", j), times = n()),
               "precision" = rep(methods[[i]], times = n()) |>
                 as_factor())
      p <- df |>   
        ggplot(aes(Var1, Var2, fill = value)) +
        geom_tile(color = "black")+
        coord_flip()+
        facet_wrap(class ~ precision, nrow = k, 
                   labeller = function (labels) {
                     labels <- lapply(labels, as.character)
                     a <-  do.call(paste, c(labels, list(sep = ",")))
                     list(gsub("\\,","-",a))
                   })+
        scale_y_continuous(breaks = scales::breaks_pretty())+
        scale_x_continuous(breaks = scales::breaks_pretty())+
        scale_fill_gradient2(low = colors[j,1], high = colors[j, 2], mid = "white", midpoint = 0, 
                             breaks = scales::breaks_pretty(n = 4)
                               # if(j == 2 & i %in% c(1, 2, 4)){
                               # scales::breaks_pretty(n = 2)
                               # } else {
                               #   scales::breaks_pretty(n = 4)
                               #   }
                             )+
        labs(fill = "", x = "Feature", y = "Feature")+
        theme(strip.background = element_rect(fill="lightgrey"), 
              legend.position = "bottom", 
              legend.box.spacing = unit(0, "pt"),# The spacing between the plotting area and the legend box (unit)
              legend.margin=margin(0,0,0,0), 
              legend.key.height = unit(.25, 'cm'), 
              legend.key.width = unit(.5, 'cm'),
              legend.text = element_text(size=7))
      
      if(i != 1) {
        p <- p + 
          theme(axis.text.y = element_blank(), 
                axis.title.y = element_blank())
      }
      if(j == 1){
        p <- p + 
          theme(axis.text.x = element_blank(), 
                axis.title.x = element_blank())
      }
      p
    }) |> 
      do.call(ggarrange_ncol1, args = _)
  })
  
  precs_plot <- ggpubr::ggarrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], plots[[6]],
                    widths = c(1.25, 1, 1, 1, 1,1), nrow = 1)
  
  vects_plot <- 
    lapply(seq_along(eig_vects), function(i){
      lapply(seq_along(eig_vects[[i]]), function(j){
        Re(eig_vects[[i]][[j]]) |> 
          # round(5) |> 
          # zapsmall() |> 
          upper_tri() |>
          reshape2::melt(na.rm = TRUE) |> 
          # filter(Var1 != Var2) |> 
          mutate("class" = rep(paste0("Class", j), times = n()),
                 "precision" = rep(methods[[i]], times = n()) |> 
                   as_factor())
      }) |> 
        do.call(rbind, args = _)
    }) |> 
    do.call(rbind, args = _) |> 
    ggplot(aes(Var1, Var2, fill = value)) +
    geom_tile(color = "black")+
    coord_flip()+
    facet_wrap(class ~ precision, nrow = k, 
               labeller = function (labels) {
                 labels <- lapply(labels, as.character)
                 a <-  do.call(paste, c(labels, list(sep = ",")))
                 list(gsub("\\,","-",a))
               })+
    scale_y_continuous(breaks = scales::breaks_pretty())+
    scale_x_continuous(breaks = scales::breaks_pretty())+
    scale_fill_gradient2(low = "red", high = "#52c569", mid = "white", midpoint = 0)+
    theme(strip.background = element_rect(fill="lightgrey"), legend.position = "bottom")+
    labs(y = "Eigen-Vector", x = "Eigen-Vector", fill = "Magnitude")
  
  vals_plot <- eig_vals |> 
    pivot_longer(starts_with("class"), names_to = "class") |> 
    mutate("class_prec" = paste0(class, "-", precision)) |> 
    ggplot(aes(Re(vector), Re(value), color = precision))+
    geom_point()+
    geom_line()+
    facet_wrap(class ~ precision, nrow = k, scales = "free_y", 
               labeller = function (labels) {
                 labels <- lapply(labels, as.character)
                 a <-  do.call(paste, c(labels, list(sep = ",")))
                 list(gsub("\\,.*"," ",a))
               })+
    theme(strip.background = element_rect(fill="lightgrey"), 
          axis.text.y = element_text(size = 7), 
          legend.position = "bottom")+
    scale_color_manual(values = c("purple", "#6495ED", "#21918c", "#5ec962", "#fde725", "black"))+
    labs(x = "Eigen-Vector", y = "Eigen-Value", color = "Precision")
  
  list(precs_plot, vects_plot, vals_plot)
}

# M Wrapper Fn ------------------------------------------------------------
sdr_mats_sim <- function(data, lambdas, nsims, dims, gamma = NULL, type = "L1", standardize_xbar = FALSE){
  x <- now()
  sim <- parallel::mclapply(1:nsims, function(sim){
    set.seed(sim)
    # 10-fold cv split
    splits <- (data |> 
                 mutate("n" = 1:n()) |> 
                 relocate(n, .after = class) |> 
                 rsample::vfold_cv(strata = class))$splits
    train <- lapply(splits, rsample::training)
    test <- lapply(splits, rsample::testing)
    
    methods <- list(NULL, "Haff", "RidgeShrinkage", "BGP16" , "SCPME")
    
    # SYS classifier model with training
    results <- lapply(seq_along(train), function(j){
      lapply(seq_along(methods), function(i){
        sdr(class~.-n, data = train[[j]], method = "SDRS", prec.est = methods[[i]], fitted.class = FALSE, 
             lam = lambdas, gamma = gamma, 
             type = type, standardize_xbar = standardize_xbar)
      })
    })
    
    Ms <- lapply(seq_along(methods), function(i){
      (lapply(seq_along(results), function(j){
        results[[j]][[i]]$M
      }) |> 
        Reduce("+", x = _))/10
    })
    
    ProjMats <- lapply(seq_along(methods), function(i){
      (lapply(seq_along(results), function(j){
        results[[j]][[i]]$ProjectionMatrix
      }) |>
        Reduce("+", x = _))/10
    })
    
    proj_dat <- lapply(seq_along(methods), function(i){
      lapply(seq_along(results), function(j){
        project(test[[j]][,-c(1:2)], results[[j]][[i]]$ProjectionMatrix[1:dims[[i]],]) |> 
          cbind("obs" = test[[j]]$n, "class" = test[[j]]$class)
      }) |> 
        do.call(rbind, args = _)
    })
    
    list(Ms, ProjMats, proj_dat)
  }, mc.cores = parallel::detectCores() - 1) 
  
  out <- lapply(1:5, function(i) {
    #Ms
    list("M" =
      (lapply(seq_along(sim), function(j) {
      sim[[j]][[1]][[i]]
    }) |>
      Reduce("+", x = _)) / nsims,
    #ProjMats
    "ProjMat" =
    (lapply(seq_along(sim), function(j) {
      sim[[j]][[2]][[i]]
    }) |>
      Reduce("+", x = _)) / nsims,
    #ProjDat
    "ProjDat" = 
    lapply(seq_along(sim), function(j){
      sim[[j]][[3]][[i]]
    }) |> 
      do.call(rbind, args = _) |> 
      as.data.frame() |> 
      group_by(obs) |> 
      group_split() |> 
      lapply(colMeans) |> 
      do.call(bind_rows, args = _) |> 
      select(-obs)
    )
  })
  
  # out <- lapply(1:5, function(i){
  #       (lapply(seq_along(sim), function(j) {
  #       sim[[j]][[i]]
  #     }) |>
  #       Reduce("+", x = _)) / nsims
  # })
  
  names(out) <- c("MLE", "Haff", "Wang", "Bodnar", "MRY")
  message(paste0("Time Difference of ", now() - x))
  system("say Just finished!")
  out
}

# Sim Wrapper Fns ---------------------------------------------------------
qda_shrink_sim <- function(data, lambdas, nsims, gamma = NULL, type = "L1", standardize_xbar = FALSE){
  x <- now()
  sim <- parallel::mclapply(1:nsims, function(sim){
    set.seed(sim)
    # 10-fold cv split
    splits <- (data |> rsample::vfold_cv(strata = class))$splits
    train <- lapply(splits, rsample::training)
    test <- lapply(splits, rsample::testing)
    
    methods <- list(NULL, "Haff", "RidgeShrinkage", "BGP16" , "SCPME")
    
    # SYS classifier model with training
    results <- lapply(seq_along(train), function(j){
      lapply(seq_along(methods), function(i){
        qda_shrink(train[[j]][,-1], train[[j]]$class, xtest = test[[j]][,-1], 
                   prec.est = methods[[i]], lam = lambdas, gamma = gamma, 
                   type = type, standardize_xbar = standardize_xbar)
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
  message(paste0("Time Difference of ", now() - x))
  ## Return
  system("say Just finished!")
  sim
}

sdr_sim <- function(data, lambdas, nsims, gamma = NULL, type = "L1", standardize_xbar = FALSE){
  x <- now()
  sim <- parallel::mclapply(1:nsims, function(sim){
    set.seed(sim)
    # 10-fold cv split
    splits <- (data |> rsample::vfold_cv(strata = class))$splits
    train <- lapply(splits, rsample::training)
    test <- lapply(splits, rsample::testing)
    
    methods <- list(NULL, "Haff", "RidgeShrinkage", "BGP16" , "SCPME")
    
    # SYS classifier model with training
    mods <- lapply(seq_along(train), function(j){
      lapply(seq_along(methods), function(i){
        sdr(class~., data = train[[j]], method = "SDRS", prec.est = methods[[i]], 
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
  system("say Just finished!")
  sim
}

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
  
  set.seed(1)
  msda_cv <- TULIP::cv.msda(as.matrix(data[,-1]), y = as.numeric(data$class), nfolds = 10)

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
  
  system("/Applications/MATLAB_R2023b.app/bin/matlab -nodisplay -batch \"run('~/Library/CloudStorage/GoogleDrive-derik.todd@gmail.com/My Drive/Research/PPP/R/prec_est/real_data_sims/ends.m'); exit\"")

  ends_error <- read_csv(file.path(sim_data_pth, "ENDS_error.csv"),
           col_names = FALSE,
           show_col_types = FALSE)[[1]]
  
  message(now() - x)
  cbind(errors, "ends" = ends_error)
}

sim_sdr_boxplot <- function(sim, dims = NULL){
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





# # extra -------------------------------------------------------------------
# 
# 
# 
# 
# y <- now()
# 
# # S <- S[[1]]
# 
# 
# out
# now() - y
# cv_qda <- function(x, grouping, xtest = NULL, lambda, A, B, C, folds = 10, nsims = NULL){
#  
#   if(is.null(folds)){
#     x <- as.matrix(x)
#     classes <- unique(grouping)
#     data <- lapply(1:length(classes), function(i) {
#       x[grouping == classes[[i]],]
#     })
#     
#     p <- ncol(x)
#     n <- lapply(data, nrow)
#     N <- length(grouping)
#     
#     priors <- lapply(n, function(x)
#       x / N)
#     S <- lapply(data, stats::cov)
#     xbar <- lapply(data, colMeans)
#     
#     B <- lapply(xbar, function(x) cbind(matrix(x, ncol = 1), 1*diag(p)))
#     A <- lapply(B, function(x) t(x))
#     C <- matrix(0, nrow = nrow(A[[1]]), ncol = ncol(B[[1]]))
#     initOmega = lapply(S,function(x) diag(diag(x)^(-1)))
#     init = lapply(seq_along(S), function(k) A[[k]] %*% initOmega[[k]] %*% B[[k]] - C)
#     zeros = matrix(0, nrow = nrow(C), ncol = ncol(C))
#     
#     S_inv <- vector("list", length(data))
#     for (k in seq_along(data)) {
#       S_inv[[k]] <- SCPME::ADMMc(S[[k]], A[[k]], B[[k]], C, initOmega[[k]], initZ = init[[k]], initY = zeros, lam = lambda[[k]])$Omega
#     }
#     
#     if (!is.null(xtest)) {
#       x <- as.matrix(xtest)
#       N <- nrow(x)
#     }
#     
#     d <- lapply(seq_along(data), function(k) {
#       # logs <- log(det(S[[i]])) - 2*log(priors[[i]])
#       out <- vector(length = N)
#       for (i in 1:N) {
#         diff <- matrix((x[i, ] - xbar[[k]]), ncol = 1)
#         out[[i]] <- -.5 * (t(diff) %*% S_inv[[k]] %*% diff - determinant(S_inv[[k]])$mod) + log(priors[[k]])
#       }
#       out
#     }) |>
#       do.call(cbind, args = _)
#     
#     classes[apply(d, 1, which.max)]
#   } else {
#     data <- cbind(grouping, x) |> as.data.frame()
#     
#     # set.seed(sim)
#     splits <- (rsample::vfold_cv(data, strata = grouping, v = 10))$splits
#     train <- lapply(splits, rsample::training)
#     test <- lapply(splits, rsample::testing)
#     
#     x_train <- lapply(train, function(x) as.matrix(x[,-1]))
#     classes <- lapply(train, function(x) unique(x$grouping))
#     data <- lapply(seq_along(train), function(i) {
#       lapply(1:length(classes[[i]]), function(k) {
#         x_train[[i]][train[[i]]$grouping == classes[[i]][[k]], ]
#       })
#     })
#       
#     
#     p <- ncol(x)
#     n <- lapply(seq_along(data), function(i){
#       lapply(data[[i]], nrow)
#     })
#                 
#     N <- lapply(seq_along(x_train), function(i){
#       nrow(x_train[[i]])
#     })
#     
#     priors <- 
#       lapply(seq_along(n), function(i){
#         lapply(n[[i]], function(x)
#           x / N[[i]])
#       })
#     S <- lapply(seq_along(x_train), function(i){
#       lapply(data[[i]], stats::cov)
#     })
#     xbar <- lapply(seq_along(x_train), function(i){
#       lapply(data[[i]], colMeans)
#     })
#     
#     # B <- lapply(seq_along(x_train), function(i){
#     #   lapply(xbar[[i]], function(x) cbind(matrix(x, ncol = 1), 1*diag(p)))
#     # })
#     # A <- lapply(seq_along(B), function(i){
#     #   lapply(B[[i]], function(x) t(x))
#     # })
#     # C <- matrix(0, nrow = nrow(A[[1]][[1]]), ncol = ncol(B[[1]][[1]]))
#     initOmega <-  lapply(seq_along(S), function(i){
#       lapply(S[[i]],function(x) diag(diag(x)^(-1)))
#     })
#     
#     B <- diag(p)
#     A <- diag(p)
#     C <- matrix(0, nrow = nrow(A), ncol = ncol(B))
#     # initZ <- lapply(seq_along(S), function(i){
#     #   lapply(seq_along(S[[i]]), function(k) A[[i]][[k]] %*% initOmega[[i]][[k]] %*% B[[i]][[k]] - C)
#     # })
#     initZ <- lapply(seq_along(S), function(i){
#       lapply(seq_along(S[[i]]), function(k) A %*% initOmega[[i]][[k]] %*% B - C)
#     })
#     zeros <-  matrix(0, nrow = nrow(C), ncol = ncol(C))
#     
#     if(is.vector(lambda)) lambda <- matrix(lambda, nrow = 1)
#     
#     parallel::mclapply(1:nrow(lambda), function(i){
#       
#       S_inv <- lapply(seq_along(data), function(v){
#         lapply(seq_along(data[[v]]), function(k){
#           SCPME::ADMMc(S[[v]][[k]], A, B, C, initOmega[[v]][[k]], initZ[[v]][[k]], zeros, lam = lambda[i, k])$Omega |> 
#             as.matrix()
#         })
#       })
#       
#       error <- lapply(seq_along(data), function(v){
#           d <- lapply(seq_along(data[[v]]), function(k) {
#             # logs <- log(det(S[[i]])) - 2*log(priors[[i]])
#             out <- vector(length = nrow(test[[v]]))
#             for (i in 1:length(out)) {
#               diff <- t(as.matrix((test[[v]][i, -1] - xbar[[v]][[k]])))
#               out[[i]] <- -.5 * (t(diff) %*% S_inv[[v]][[k]] %*% diff - determinant(S_inv[[v]][[k]])$mod) + log(priors[[v]][[k]])
#             }
#             out
#           }) |>
#             do.call(cbind, args = _)
#           mean(test[[v]][,1] != classes[[v]][apply(d, 1, which.max)])
#         }) |> 
#           do.call(rbind, args = _) |> 
#           colMeans()
#       matrix(c(lambda[i, ], error), nrow = 1)
#     }, mc.cores = parallel::detectCores() - 1) |> 
#       do.call(rbind, args = _)
#   }  
# }
# 
# x <- now()
# y <- cv_qda(dat[,-1], dat$class, lambda = lams)
# now() - x
# 
# 
# z[1:10,3] == (y[1:10,3] |> unlist())
