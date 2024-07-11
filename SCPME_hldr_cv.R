# qda ---------------------------------------------------------------------
SCPME_qda_cv <- function(lambdas, data, nsims = 1, gamma = NULL){
  
  if(is.vector(lambdas)) {nlams <- 1; lambdas <- matrix(lambdas, nrow = 1)}
  if(!is.vector(lambdas)) nlams <- nrow(lambdas)
  
  parallel::mclapply(1:nlams, function(k){
    error <-  lapply(1:nsims, function(sim){
      set.seed(sim)
      splits <- (data |> rsample::vfold_cv(strata = class))$splits
      train <- lapply(splits, rsample::training)
      test <- lapply(splits, rsample::testing)
      
      
      # SYS classifier model with training
      results <- lapply(seq_along(train), function(j){
        qda_shrink(train[[j]][,-1], train[[j]]$class, xtest = test[[j]][,-1], prec.est = "SCPME", lam = lambdas[k,], gamma = gamma, 
                   type = "qda",
                   standardize_xbar = TRUE
                   )
      })
      
    lapply(seq_along(results), function(i){
          mean(results[[i]] != test[[i]]$class)
        }) |> 
          unlist() |> 
          mean()
    }) 
      
    meds <- median(error[[1]][[1]])
    out <- matrix(c(lambdas[k,], meds, error[[1]][[2]]), nrow = 1)
    colnames(out) <- c(paste0("lambda", 1:ncol(lambdas)), "error", "lik")
    if(nsims > 1) {
      out <- cbind(out, "sd" = sd(error))
    }
    out
    
  }, mc.cores = parallel::detectCores() - 1) |> 
    do.call(rbind, args = _)
}


# hldr --------------------------------------------------------------------
SCPME_hldr_cv <- function(lambdas, data, nsims = 1, gamma = NULL){
  if(is.vector(lambdas)) {nlams <- 1; lambdas <- matrix(lambdas, nrow = 1)}
  if(!is.vector(lambdas)) nlams <- nrow(lambdas)
  
  parallel::mclapply(1:nlams, function(k){
   error <-  lapply(1:nsims, function(sim){
      set.seed(sim)
      splits <- (dat |> rsample::vfold_cv(strata = class))$splits
      train <- lapply(splits, rsample::training)
      test <- lapply(splits, rsample::testing)
      
      methods <- list(NULL, "Haff", "BGP16", "RidgeShrinkage", "SCPME")
      
      # SYS classifier model with training
      mods <- lapply(seq_along(train), function(j){
        hldr(class~., data = train[[j]], method = "SYS", prec.est = methods[[5]], lam = unlist(lambdas[k,]), fitted.class = FALSE, gamma = gamma, 
             type = "qda",
             standardize_xbar = TRUE
             )
      })
      
      results <- lapply(seq_along(mods), function(j){
        lapply(1:(ncol(train[[j]]) - 2), function(r){
          predict(mods[[j]], dims = 1:r,
                  newdata = test[[j]][,-1], grouping = test[[j]]$class)$class
        })
      })
      
      lapply(seq_along(mods), function(i){
        lapply(seq_along(results[[i]]), function(j){
          mean(results[[i]][[j]] != test[[i]]$class)
        }) |> 
          do.call(cbind, args = _)
      }) |> 
        do.call(rbind, args = _) |> 
        colMeans()
    }) |> 
      do.call(rbind, args = _)
   meds <- apply(error, 2, median)
   min_dims <- which.min(meds)
   out <- matrix(c(lambdas[k,], min_dims, meds[min_dims]), nrow = 1)
   colnames(out) <- c(paste0("lambda", 1:ncol(lambdas)), "ndims", "error")
   if(nsims > 1) {
     out <- cbind(out, "sd" = apply(error, 2, sd)[min_dims])
   }
   out
  }, mc.cores = parallel::detectCores() - 1) |> 
    do.call(rbind, args = _)
}


# -------------------------------------------------------------------------
SCPME_hldr_pca_cv <- function(x, grouping, lambdas, npca, dimselect_cv = NULL, nsims = 1, gamma = NULL){
  if(is.vector(lambdas)) {nlams <- 1; lambdas <- matrix(lambdas, nrow = 1)}
  if(!is.vector(lambdas)) nlams <- nrow(lambdas)
  
  
  
  parallel::mclapply(1:nlams, function(k){
    
  if(!is.null(dimselect_cv)){
    crit_cv <- dimselect_cv(x = x, grouping = grouping, cv_method  = dimselect_cv, npca = npca, nsims = 10, lambdas = unlist(lambdas[k,]))
    dims_cv <- sort(crit_cv, decreasing = ifelse(dimselect_cv == "t", TRUE, FALSE), index.return = TRUE)$ix
  } else{
    dims_cv <- 1:npca
  }
    
    error <-  lapply(1:nsims, function(sim){
      set.seed(sim)
      splits <- (dat |> rsample::vfold_cv(strata = class))$splits
      train <- lapply(splits, rsample::training)
      test <- lapply(splits, rsample::testing)
      
      lapply(seq_along(train), function(j){
        ## prin comp projection
        train <- train[[j]]
        test <- test[[j]]
        # Perform PCA on the training data
        pca_result <- prcomp(train[,-1], center = TRUE, scale. = TRUE)
        
        # pca_dims_ord <- dimselect(train[,-1], train$class, ProjectionMatrix = t(pca_result$rotation), folds = 10, method = "t")$dims
        
        # (pca_result$sdev[pca_dims_ord][1:36]^2/sum((pca_result$sdev^2)[1:36]) |> 
        #   cumsum())
        # npca <- -(summary(pca_result)$importance[2,][pca_dims_ord] |> 
        #             cumsum())[1:36] |> 
        #   lik() |> 
        #   which.max()
        # Print PCA summary
        # summary(pca_result)
        
        
        # Print PCA loadings (projection matrix)
        # pca_result$rotation
        # The projection matrix is the rotation matrix
        pca_proj <- pca_result$rotation[,1:npca]
        # Standardize the test data using the training data parameters
        # Extract the mean and standard deviation from the training data scaling
        train_mean <- colMeans(train[,-1])
        train_sd <- apply(train[,-1], 2, sd)
        
        # Apply the same scaling to the test data
        scaled_test_data <- scale(test[,-1], center = train_mean, scale = train_sd)
        # Choose the number of principal components to use (e.g., the first 2)
        # num_components <- 36
        df <- bind_cols("class" = train$class, as.matrix(train[,-1]) %*% pca_proj)
        
        
        mod_hldr <- hldr::hldr(class ~., data = df, method = "SYS", prec.est = "SCPME", lam = unlist(lambdas[k,]), dims = dims_cv)
      
      
      ## cv to select dims by T
        
        
      lapply(1:npca, function(r){
          pred <- predict(mod_hldr, newdata = as.matrix(scaled_test_data) %*% pca_proj, dims = 1:r)$class
          mean(pred != test$class)
        }) |> 
        do.call(cbind, args = _)
      }) |> 
        do.call(rbind, args = _) |> 
        colMeans()
    }) |> 
      do.call(rbind, args = _)
    meds <- apply(error, 2, median)
    min_dims <- which.min(meds)
    out <- matrix(c(lambdas[k,], min_dims, meds[min_dims]), nrow = 1)
    colnames(out) <- c(paste0("lambda", 1:ncol(lambdas)), "ndims", "error")
    if(nsims > 1) {
      out <- cbind(out, "sd" = apply(error, 2, sd)[min_dims])
    }
    out
  }, mc.cores = parallel::detectCores() - 1) |> 
    do.call(rbind, args = _)
}
# Examples ----------------------------------------------------------------