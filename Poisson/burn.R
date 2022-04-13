#GH(b, Sigmai, gh)

gh <- 3
GH <- statmod::gauss.quad(gh, 'hermite')
b.hat <- b.hat; n <- length(b.hat)
q <- 2
C <- lapply(Sigmai, chol)
# b.var <- mapply(function(Z, S) tcrossprod(Z %*% S, Z), S = Sigmai, Z = Z)
b.hat.new <- vector('list', n)
v <- as.matrix(expand.grid(lapply(1:q, function(k, u) u$nodes, u = GH)))
for(i in 1:n){
  b.hat.new[[i]] <- t(sqrt(2) * solve(C[[i]], t(v)) + b.hat[[i]])
}
w <- as.matrix(expand.grid(lapply(1:q, function(k, u) u$weights, u = GH)))
w <- 2^(q/2) * apply(w,1,prod) * exp(rowSums(v * v))
b2 <- lapply(b.hat.new, function(b) t(apply(b,1,tcrossprod)))
Zb <- mapply(function(Z, b) exp(tcrossprod(Z, b)) %*% w, Z = Z, b = b.hat.new)

eta.quad <- mapply(function(X, Z) X %*% beta + Z, X = X, Z = Zb)

eta <- mapply(function(X, Z, b) X %*% beta + Z %*% b, X = X, Z = Z, b = b.hat, SIMPLIFY = F)
sigma <- mapply(function(Z, S) unname(sqrt(diag(tcrossprod(Z %*% S, Z)))), Z = Z, S = Sigmai)

statmod::gauss.quad.prob(3, 'normal', mu = eta[[1]], sigma = sigma[[1]])


GHfun <- function (b, y_lis, N_lis, X_lis, Z_lis, offset_lis, X_zi_lis, Z_zi_lis, offset_zi_lis,
                   betas, inv_D, phis, gammas, k, q,
                   canonical, user_defined, Zty_lis, log_dens, mu_fun, var_fun, mu.eta_fun,
                   score_eta_fun, score_phis_fun, score_eta_zi_fun) {
  GH <- gauher(k)
  aGH <- find_modes(b, y_lis, N_lis, X_lis, Z_lis, offset_lis, X_zi_lis, Z_zi_lis, 
                    offset_zi_lis, betas, inv_D, phis, gammas,
                    canonical, user_defined, Zty_lis, log_dens, mu_fun, var_fun, 
                    mu.eta_fun, score_eta_fun, score_phis_fun, score_eta_zi_fun)
  modes <- aGH$post_modes
  chol_hessians <- lapply(aGH$post_hessian, chol)
  b <- as.matrix(expand.grid(lapply(seq_len(q), function (k, u) u$x, u = GH)))
  n <- nrow(modes)
  b_new <- vector("list", n)
  log_dets <- numeric(n)
  for (i in seq_len(n)) {
    b_new[[i]] <- t(sqrt(2) * solve(chol_hessians[[i]], t(b)) + modes[i, ])
    log_dets[i] <- - determinant.matrix(chol_hessians[[i]], logarithm = TRUE)$modulus
  }
  wGH <- as.matrix(expand.grid(lapply(seq_len(q), function (k, u) u$w, u = GH)))
  wGH <- 2^(q/2) * apply(wGH, 1, prod) * exp(rowSums(b * b))
  b2 <- lapply(b_new, function (b) if (q == 1) b * b else
    t(apply(b, 1, function (x) x %o% x)))
  ind_Z <- seq_len(ncol(Z_lis[[1]]))
  Ztb <- do.call('rbind', mapply(function (z, b) z %*% t(b[, ind_Z, drop = FALSE]), 
                                 Z_lis, b_new, SIMPLIFY = FALSE))
  Z_zitb <- if (!is.null(Z_zi_lis[[1]])) {
    do.call('rbind', mapply(function (z, b) z %*% t(b[, -ind_Z, drop = FALSE]), 
                            Z_zi_lis, b_new, SIMPLIFY = FALSE))  
  } 
  list(b = do.call('rbind', b_new), b2 = do.call('rbind', b2), Ztb = Ztb, Z_zitb = Z_zitb,
       wGH = wGH, log_dets = log_dets, post_modes = modes, 
       post_vars = lapply(aGH$post_hessian, solve))
}