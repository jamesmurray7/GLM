# Investigating why certain Hdelta calls fails (choice of epsilon?)

Sdelta$Hessian[to.repeat]
start.eps <- .Machine$double.eps^.25

sink('~/Desktop/hess-test.txt')
invisible(sapply(to.repeat, function(i){
  cat(sprintf("i: %d, Hessian on pow(datum::eps, 1./4.): %.3f\n", i, Sdelta$Hessian[i]))
  # A) Smaller eps?
  H.new1 <- Hdelta(delta[[i]], b[[i]], X[[i]], Z[[i]], Y[[i]],
                  lY[[i]], beta, tau[[i]], w, v, summax[[i]], .Machine$double.eps^(1/3))
  cat(sprintf("... on pow(datum::eps, 1./3.): %.3f\n", H.new1))
  
  # B) Even smaller!
  H.new2 <- Hdelta(delta[[i]], b[[i]], X[[i]], Z[[i]], Y[[i]],
                   lY[[i]], beta, tau[[i]], w, v, summax[[i]], .Machine$double.eps^(1/2))
  cat(sprintf("... on pow(datum::eps, 1./2.): %.3f\n", H.new2))
  
  
  # C) Bigger eps
  H.new3 <- Hdelta(delta[[i]], b[[i]], X[[i]], Z[[i]], Y[[i]],
                  lY[[i]], beta, tau[[i]], w, v, summax[[i]], 1e-3)
  cat(sprintf("... on 1e-3: %.3f\n", H.new3))
  
  # D) even bigger...
  H.new4 <- Hdelta(delta[[i]], b[[i]], X[[i]], Z[[i]], Y[[i]],
                   lY[[i]], beta, tau[[i]], w, v, summax[[i]], 1e-2)
  cat(sprintf("... on 1e-2: %.3f\n", H.new4))
  
  
  # Were either good?
  if(H.new1 < 0) cat('Smaller eps worked!\n')
  if(H.new2 < 0) cat('Smallest eps worked!\n')
  if(H.new3 < 0) cat('Larger eps worked!\n')
  if(H.new4 < 0) cat('largest eps worked!\n')
  if(H.new1 > 0 && H.new2 > 0 && H.new3 > 0 && H.new4 > 0) cat('Nothing Worked!\n')
  
  cat('\n\n')
  return(invisible(1))
}))
sink()

# Trying another way ...


eta <- X[[1]] %*% beta + Z[[1]] %*% b.hat[[1]]
mu <- exp(eta)
nu <- rep(exp(delta[[1]]), length(mu))
xi <- summax[[1]]
loglam <- log(lambda_appx(mu, nu, xi))
logZ <- logZ_c(loglam, nu, xi)

# Expected value and variances.
YlY <- E_YlY(loglam, logZ, nu, xi)
lY <- E_lY(loglam, logZ, nu, xi)
VY <- V_Y(loglam, logZ, nu, xi)
VlY <- V_lY(loglam, logZ, nu, lY, xi)
y <- Y[[1]]

A <- YlY - mu * lY
S <- sum((A * (y - mu) / VY - lgamma(y+1) + lY) * nu)
H <- -sum(((-(A)^2 / VY + VlY) * nu^2))
delta[[1]]-S/H

#' At point estimate for \hat{\b} 
sapply(inds.met, function(i){
  y <- Y[[i]]
  # Means, rates, dispersion etc
  eta <- X[[i]] %*% beta + Z[[i]] %*% b.hat[[i]]
  mu <- exp(eta)
  nu <- rep(exp(delta[[i]]), length(mu))
  xi <- summax[[i]]
  loglam <- log(lambda_appx(mu, nu, xi))
  logZ <- logZ_c(loglam, nu, xi)
  # Expected value and variances.
  YlY <- E_YlY(loglam, logZ, nu, xi)
  lY <- E_lY(loglam, logZ, nu, xi)
  VY <- V_Y(loglam, logZ, nu, xi)
  VlY <- V_lY(loglam, logZ, nu, lY, xi)
  A <- YlY - mu * lY
  S <- sum((A * (y - mu) / VY - lgamma(y+1) + lY) * nu)
  H <- -sum(((-(A)^2 / VY + VlY) * nu^2))
  delta[[i]]-S/H
}) -> new


#' With quadrature
new.q <- sapply(inds.met, function(i){
  y <- Y[[i]]
  # Linear predictor
  eta <- X[[i]] %*% beta + Z[[i]] %*% b.hat[[i]]
  nu <- rep(exp(delta[[i]]), length(eta))
  xi <- summax[[i]]
  S <- H <- 0
  for(l in 1:3){
    eta.l <- eta + tau[[i]] * v[l]
    mu <- exp(eta.l)
    loglam <- log(lambda_appx(mu, nu, xi))
    logZ <- logZ_c(loglam, nu, xi)
    # Expected value and variances.
    YlY <- E_YlY(loglam, logZ, nu, xi)
    lY <- E_lY(loglam, logZ, nu, xi)
    VY <- V_Y(loglam, logZ, nu, xi)
    VlY <- V_lY(loglam, logZ, nu, lY, xi)
    A <- YlY - mu * lY
    S <- S + w[l] * sum((A * (y - mu) / VY - lgamma(y+1) + lY) * nu)
    H <- H + w[l] * -sum(((-(A)^2 / VY + VlY) * nu^2))
  }
  delta[[i]] - S/H
})

new.cpp <- sapply(inds.met, function(i){
  new_delta_update(delta[[i]], X[[i]], Z[[i]], Y[[i]],
                   b.hat[[i]], beta, summax[[i]], w, v, tau[[i]], F)
})
new.q.cpp <- sapply(inds.met, function(i){
  new_delta_update(delta[[i]], X[[i]], Z[[i]], Y[[i]],
                   b.hat[[i]], beta, summax[[i]], w, v, tau[[i]], T)
})
  