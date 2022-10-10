#' Testing possible way to obtain subject-specific SEs (and how similar they are to H(delta))...
news <- vector('list')
for(i in inds.met){
  delta.i <- delta[[i]]
  xi.i <- summax[[i]]
  lb <- if(delta.i > 0) 0 else delta.i - .5
  ub <- if(delta.i < 0) 0 else delta.i + .5
  cat(sprintf("delta: %.2f, [%.2f, %.2f]\n", delta.i, lb, ub))
  
  mu.i <- exp(X[[i]] %*% Omega$beta + Z[[i]] %*% b[[i]])
  
  optim(delta.i, ff3, NULL,
        Y = Y[[i]], mu = mu.i, summax = xi.i,
        method = 'Brent',
        lower = lb, upper = ub, hessian = T) -> o
  news[[i]] <- o
}

which.max(sapply(inds.met, function(i) news[[i]]$par - delta[[i]]))

SEs <- sapply(inds.met, function(x) sqrt(solve(news[[x]]$hessian)))
delts <- sapply(inds.met, function(x) news[[x]]$par)

delta.df2 <- data.frame(
  id = inds.met, 
  delta = delts,
  lb = delts - 1.96 * SEs,
  ub = delts + 1.96 * SEs
)
