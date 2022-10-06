mus <- mapply(function(x,z,b) exp(x %*% beta + z %*% b), x = X, z = Z, b = b, SIMPLIFY = F)
nus <- mapply(function(d, y) rep(exp(d), length(y)), d = delta, y = Y, SIMPLIFY = F)

lam <- lambda_appx(mus[[1]], nus[[1]], summax[[1]])
loglam <- log(lam)
logZ <- logZ_c(loglam, nus[[1]], summax[[1]])



out <-  matrix(0, nr = length(nus[[1]]), nc = summax[[1]] + 1)
for(i in 0:summax[[1]]){
  out[,(i+1)] <- (i - mus[[1]])^2 * lam^i / (factorial(i)^nus[[1]] * exp(logZ))
}
out
apply(out, 1, sum)
V_Y(loglam, logZ, nus[[1]], summax[[1]])


mu <- exp(eta)
loglam <- log(lambda_appx(mu, nu, xi))
logZ <- logZ_c(loglam, nu, xi)
# Expected value and variances.
YlY <- E_YlY(loglam, logZ, nu, xi)
lY <- E_lY(loglam, logZ, nu, xi)
VY <- V_Y(loglam, logZ, nu, xi)
VlY <- V_lY(loglam, logZ, nu, lY, xi)
A <- YlY - mu * lY

Snu <- A * (Y[[1]] - mu) / VY - lgamma(Y[[1]] + 1) + lY

S <- sum((A * (Y[[1]] - mu) / VY - lgamma(Y[[1]]+1) + lY) * nu)
rhs <- diag(nu) %*% matrix(1, nr = 9, nc= 1)
crossprod(Snu, rhs)

H <- -sum(((-(A)^2 / VY + VlY) * nu^2))

W <- matrix(1, nr = 9, nc= 1)
diag(c(diag(nu) %*% W)) %*% W
crossprod(-(A)^2 / VY + VlY, diag(c(diag(nu^2) %*% W)) %*% W)




# Testing Huang's variance vs Fung/Shmueli --------------------------------
# Huang's includes mean, rate, disp
# Shmueli includes only rate, disp.

fn <- function(mus, nus, summax){
  lam <- lambda_appx(mus, nus, summax)
  loglam <- log(lam)
  logZ <- logZ_c(loglam, nus, summax)
  
  
  
  out <-  matrix(0, nr = length(nus), nc = summax + 1)
  for(i in 0:summax){
    out[,(i+1)] <- (i - mus)^2 * lam^i / (factorial(i)^nus * exp(logZ))
  }
  out
  a <- apply(out, 1, sum)
  b <- V_Y(loglam, logZ, nus, summax)
  list(a=a,b=b)
}

test <- mapply(function(m,n,x){
  fn(m,n,x)
}, m = mus, n = nus, x = summax, SIMPLIFY = F)
y <- Y[[1]]
tests <- vector('list', 2*max(y)-max(y))
for(s in seq_along(seq(max(y), 2 * max(y)))){
  tests[[s]] <- fn(mus[[1]], nus[[1]], seq(max(y), 2 * max(y))[s])
}

ylim <- c(min(sapply(tests, function(x) min(x$b))),
          max(sapply(tests, function(x) max(x$b))))
xlim <- c(min(sapply(tests, function(x) min(x$a))),
          max(sapply(tests, function(x) max(x$a))))
plot(tests[[13]]$a, tests[[13]]$b, pch = 20,
     ylim = ylim, xlim = xlim, xlab = 'Huang', ylab = 'Fung/Shmueli')
abline(0,1)
for(s in seq_along(seq(max(y) + 1, 2 * max(y)))){
  points(tests[[s]]$a, tests[[s]]$b, col = s, pch = 20)
}
legend('topleft',
       pch = 20, ncol = 3,
       col = 1:length(seq(max(y) + 1, 2 * max(y))),
       legend = as.character( seq(max(y) + 1, 2 * max(y))),
       bty = 'n', cex = .5) # Very little difference after xi = 2*max(Y).

       