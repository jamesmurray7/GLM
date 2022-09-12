ff <- function(x){
  joint_density(c(x), X[[1]], Y[[1]], lY[[1]], Z[[1]], G[[1]], beta, delta, D, S[[1]], SS[[1]],
                Fi[[1]], Fu[[1]], l0i[[1]], l0u[[1]], Delta[[1]], gamma, zeta, lambda.mat, V.mat, logZ.mat,
                all.mus, round(nu.vec, 2), summax = summax)
}
xs <- seq(-2, 2, .001)
fxs <- sapply(xs, ff)
plot(fxs ~ xs, type = 'l')
abline(v = xs[which.min(fxs)], col='magenta')
abline(v = nloptr::bobyqa(b[[1]], ff)$par, col = 'blue')
abline(v = optim(b[[1]], ff, NULL, method = 'BFGS')$par, col= 'red')
abline(v = nloptr::bobyqa(c(0), ff)$par, col = 'green')
abline(v = optim(c(0), ff, NULL, method = 'BFGS')$par, col= 'black')

b[[1]]
t1 <- nloptr::bobyqa(b[[1]], ff)
t2 <- optim(b[[1]], ff, NULL, method = 'BFGS', hessian = T)
t3 <- nloptr::bobyqa(c(0), ff)
t4 <- optim(c(0), ff, NULL, method = 'BFGS', hessian= T)
