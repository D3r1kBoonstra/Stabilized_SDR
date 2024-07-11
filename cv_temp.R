group_dat <- dat |> 
  group_by(class) |> 
  group_split() |> 
  lapply(function(x) x[,-1])

xbars <- group_dat |> 
  lapply(function(x) {
    as.matrix(colMeans(x))
  })


S_inv <- lapply(group_dat, function(x) solve(cov(x)))

lapply(xbars, mean)

eig_decomp <- lapply(S_inv, eigen)

eig_values <- lapply(eig_decomp, function(x) x$values |> as.matrix() |> zapsmall()) 

out <- SCPME_qda(dat[,-1], grouping = dat$class, lam = c(eig_values[[1]][[1]], eig_values[[2]][[1]]))
out[[1]]$Omega |> zapsmall()

# |> 
#   do.call(rbind, args = _) |> 
#   as.data.frame() |> 
#   mutate("class" = rep(1:2, each = n()/2), 
#          "vector" = rep(1:(n()/2), times = 2)) |> 
#   ggplot(aes(vector, V1, color = as_factor(class)))+
#   geom_point()+
#   geom_line()+
#   facet_wrap(~class, scales = "free_y")

smallest_eig_vec <- lapply(seq_along(eig_decomp), function(i){
  x <- eig_decomp[[i]]
  x$vectors[,3:ncol(x$vectors)] |> 
    as.matrix()
})

lm(smallest_eig_vec[[2]] ~ xbars[[2]] - 1) |> 
  summary()

coefs <- qr.solve(smallest_eig_vec[[2]], xbars[[2]])

(all.equal(smallest_eig_vec[[2]] %*% coefs, xbars[[2]]))

qr(smallest_eig_vec[[2]])$rank
qr(cbind(smallest_eig_vec[[2]], xbars[[2]]))$rank

standardize <- function(x) {(x - mean(x))/sd(x)}

B <- lapply(xbars, function(x) cbind(x, diag(nrow(x))))

lapply(xbars, function(x) mean(abs(x)))

A <- lapply(B, t)

off_diagonals <- function(matrix) {
  indices <- which(row(matrix) != col(matrix), arr.ind = TRUE)
  matrix[indices]
}


lapply(seq_along(S_inv), function(i){
  sum(abs(off_diagonals(A[[i]] %*% S_inv[[i]] %*% B[[i]])))
})

sum(diag(solve(S_inv[[1]]) %*% out[[1]]$Omega))

lapply(S_inv, function(x) sum(diag(solve(x) %*% x)))

lapply(S_inv, function(x) log(det(x)))

lams <- expand.grid(seq(.1, 2.5, by = .1), seq(.1, 4, by = .1))

x <- now()
out <- SCPME_hldr_cv(lambdas = c(1.5, 1.5), gamma = c(1.5, 3), data = dat, nsims = 10)
now() - x

glasso::glasso(solve(S_inv[[1]]), rho = eig_values[[1]][[1]])

library(plotly)
out |> 
  as_tibble() |> 
  lapply(unlist) |> 
  do.call(bind_cols, args = _) |> 
  arrange(error) |> 
  filter(error == min(error)) |>
  # arrange(sd) |> 
  # pull(lambda2)
  # filter(sd == min(sd)) |> 
  # summary()
  # filter(sd == min(sd))
  # plot_ly(x = ~lambda1, y = ~lambda2, z = ~lambda3,
  #         marker = list(color = ~error,
  #                       colorscale= "Viridis",
  #                       reversescale = TRUE,
  #                       showscale = TRUE)) |>
  # add_markers()
  ggplot(aes(lambda1, lambda2, color = sd))+
  geom_point()+
  viridis::scale_color_viridis(direction = -1)

