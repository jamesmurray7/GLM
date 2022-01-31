rm(list=ls())
source('ll.R')

vech <- function(x) x[lower.tri(x, T)]

# Data and data matrices
test <- simData_joint()
data <- test$data
n <- length(unique(data$id))
X <- Z <- G <- Y <- list()
for(i in 1:n){
  X[[i]] <- model.matrix(~ time + cont + bin, data[data$id == i, ])
  Z[[i]] <- G[[i]] <- model.matrix(~ time, data[data$id == i, ])
  Y[[i]] <- data[data$id == i, 'Y']
}

# Initial conditions
fit <- glmmTMB::glmmTMB(Y ~ time + cont + bin + (1 + time|id), 
                        dispformula = ~ time,
                        data, 
                        family = glmmTMB::nbinom2)
beta <- glmmTMB::fixef(fit)$cond
delta <- glmmTMB::fixef(fit)$disp # delta <- c(1, 0)
b <- glmmTMB::ranef(fit)$cond$id
bl <- lapply(1:n, function(i) as.numeric(b[i,]))
D <- matrix(glmmTMB::VarCorr(fit)$cond$id,2,2)
params <- c(vech(D), beta, delta)
names(params) <- c(paste0('D', apply(which(lower.tri(D,T), arr.ind = T), 1, paste0, collapse = ',')),
                   names(beta), paste0('disp_', names(delta)))


tol <- 1e-1; diff <- 100
.ngh <- 3
gh <- statmod:::gauss.quad.prob(.ngh, 'normal')
w <- gh$w;v <- gh$n
while(diff > tol) {
  #' ####
  #' E-step
  #' ####
  .bfits <- mapply(function(b, X, Y, Z, G){
    u <- ucminf::ucminf(b, .ll,  Score_b,  
                        X = X, Y = Y, Z = Z, G = G, beta = beta, delta = delta, D = D, summax = 10)
    list(u$par, u$invhessian.lt)
  }, b = bl, X = X, Y = Y, Z = Z, G = G, SIMPLIFY = F)
  message('bfits done')
  
  # Posterior mode of b ~ N(bhat, Sigma)
  b.hat <- lapply(.bfits, el, 1)
  # And posterior variance 
   Sigmai <- mapply(function(b, X, Y, Z, G){
     solve(pracma::hessian(.ll, b, X = X, Y = Y, Z = Z, G = G, beta = beta, delta = delta, D = D, summax = 10))
   }, b = b.hat, X = X, Y = Y, Z = Z, G = G, SIMPLIFY = F)
  
  Sigmai <- lapply(1:n, function(i){
    x <- Sigmai[[i]]
    if(det(x) <= 0 | any(eigen(x)$val < 0)) x <- tryCatch(as.matrix(Matrix::nearPD(x)$mat), error = function(e) NA)
    if(is.na(x)) x <- vech2mat(.bfits[[i]][[2]], 2)
    x
  })
  
  # objects for mu and nu
  mu <- mapply(function(X, Z, b){
    exp(X %*% beta + Z %*% b)
  }, X = X, Z = Z, b = b.hat, SIMPLIFY = F)
  
  nu <- mapply(function(G){
    exp(G %*% delta)
  }, G = G, SIMPLIFY = F)
  
  # lambda, from uniroot on mu and nu
  lambda <- mapply(getlambda, mu, nu, summax = 10, SIMPLIFY = F)
  log.lambda <- lapply(lambda, log)
  # Z, normalising constant using lambda and nu
  logZ_ <- mapply(logZ_c, log.lambda, nu, summax = 10, SIMPLIFY = F)
  Z_ <- lapply(logZ_, exp)
  # Means
  lfactY <- mapply(E_logfactY, lambda, nu, logZ_, summax = 10, SIMPLIFY = F)
  YlfactY <- mapply(E_YlogfactY, lambda, nu, logZ_, summax = 10, SIMPLIFY = F)
  means <- mapply(E_means, lambda, nu, logZ_, summax = 10, SIMPLIFY = F)
  V <- mapply(V_mu_lambda, mu, lambda, nu, summax = 10, SIMPLIFY = F)
  
  #' D update ----
  Drhs <- mapply(function(b, S){
    S + tcrossprod(b)
  }, b = b.hat, S = Sigmai, SIMPLIFY = F)
  
  #' beta update ----
  
  tau <- mapply(function(Z, S) unname(sqrt(diag(tcrossprod(Z %*% S, Z)))), Z = Z, S = Sigmai, SIMPLIFY = F)
  
  # Sb <- mapply(Score_beta_quad, beta = replicate(n, beta, simplify = F),
  #              Y = Y, X = X, Z = Z, b = b.hat, nu = nu, tau = tau, SIMPLIFY = F) 
  # 
  # Hb <- mapply(function(Y, X, Z, b, nu, tau){
  #   GLMMadaptive:::fd_vec(beta, Score_beta_quad, Y = Y, X = X, Z = Z, b = b, nu = nu, tau = tau)
  # }, Y = Y, X = X, Z = Z, b = b.hat, nu = nu, tau = tau, SIMPLIFY = F)
  # 
  # Hb <- Reduce('+', Hb)
  

  #' delta update ----
  # A <- mapply(function(lY, mu, YlY) YlY - mu * lY, lY = lfactY, mu = mu, YlY = YlfactY, SIMPLIFY = F)
  # B <- lfactY
  # C <- mapply(Var_logfactY, lambda, nu, logZ_, summax = 10, SIMPLIFY = F)
  
  Score_nu <- function(A, Y, mu, V, B) A * (Y - mu) / V - (lgamma(Y + 1) - B)
  
  Sd <- mapply(Score_delta, replicate(n, delta, F),
               Y = Y, X = X, Z = Z, b = b.hat, G = G, tau = tau)
  
  Hd <- mapply(function(Y, X, Z, b, G, tau){
    GLMMadaptive:::fd_vec(
      delta, Score_delta, Y = Y, X = X, Z = Z, b = b, G = G, tau = tau
    )
  }, Y = Y, X = X, Z = Z, b = b.hat, G = G, tau = tau, SIMPLIFY = F)
  
  # updates
  D.new <- D#Reduce('+', Drhs)/n
  # beta.new <- beta-solve(Hb, colSums(do.call(rbind, Sb)))
  delta.new <- delta - solve(Reduce('+', Hd), c(rowSums(Sd)))
  
  # params
  beta.new <- beta
  params.new <- c(vech(D.new), beta.new, delta.new); names(params.new) <- names(params)
  
  # Check rel. diff
  diffs <- abs(params.new-params)/(abs(params) + 1e-3)
  diff <- max(diffs)
  print(sapply(params.new, round, 3))
  message('Max relative difference: ', round(diff, 4))
  message('For ', names(params.new)[which(diffs == diff)])
  message('Old: ', round(params[which(diffs == diff)], 4), ', new -> ', round(params.new[which(diffs == diff)], 4))
  
  # Set new as old
  bl <- b.hat
  D <- D.new
  beta <- beta#.new
  delta <- delta.new
  params <- params.new
  
}
