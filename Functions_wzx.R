dahl <- function(samples) {
  n <- length(samples[[1]]$index_cluster)
  niter <- length(samples)
  z_sample <- matrix(0, nrow = niter, ncol = n)
  for (i in 1:niter) {
    z_sample[i, ] <- samples[[i]]$index_cluster
  }
  BList <- map(1:niter, ~outer(z_sample[.x,], z_sample[.x,], "=="))
  BBar <- Reduce("+", BList) / niter
  SSE <- map_dbl(BList, ~sum((.x - BBar)^2))
  return(list(min.sse = min(SSE), cLS = which.min(SSE),
              cluster = as.numeric(z_sample[which.min(SSE),])))
}


DIC <- function(samples, X_grid, B_grid, Y_grid, Dahl_index) {
  n <- length(samples[[1]]$index_cluster)
  p <- length(samples[[1]]$beta)
  q <- length(samples[[1]]$ita)
  beta_sample <- matrix(0, nrow = length(samples), ncol = p)
  ita_sample <- matrix(0, nrow = length(samples), ncol = q)
  lambda_sample <- matrix(0, nrow = length(samples), ncol = n)
  for (i in 1:length(samples)) {
    beta_sample[i, ] <- samples[[i]]$beta
    ita_sample[i, ] <- samples[[i]]$ita
    lambda_sample[i, ] <- samples[[i]]$lambda0[samples[[i]]$index_cluster]
  }
  beta_fitted <- colMeans(beta_sample)
  ita_fitted <- colMeans(ita_sample)
  lambda_fitted <- samples[[Dahl_index]]$lambda0[samples[[Dahl_index]]$index_cluster]
  ## for each poster sample draw
  dev_sample <- rep(0, length(samples))
  for (i in 1:length(samples)) {
    dev_sample[i] <- -2 * (sum(Y_grid * log(N_grid * lambda_sample[i, ] * exp(X_grid %*% beta_sample[i, ] +
                                                                                B_grid %*% ita_sample[i, ]))-
                                 N_grid * lambda_sample[i, ] * exp(X_grid %*% beta_sample[i, ] + B_grid %*% ita_sample[i, ])))
  }
  dev_estimate <- -2 * (sum(Y_grid * log(N_grid * lambda_fitted * exp(X_grid %*% beta_fitted + B_grid %*% ita_fitted))-
                              N_grid * lambda_fitted * exp(X_grid %*% beta_fitted + B_grid %*% ita_fitted)))
  pd <- mean(dev_sample) - dev_estimate
  return(list(DIC = dev_estimate + 2*pd, pD = pd))
}
