#' ####
#' sim-with-poisson
#' Simulating N lots of data which are truly over/underdispersed and fitting 
#'   using Multi-test with poisson family (i.e. ignoring dispersion...)

rm(list=ls())
setwd('../mpcmp/half_and_half/')
source('EM.R')
N <- 50
tests <- replicate(50, 
                   simData_joint2(n = 250, delta = c(1, -1), 
                        ntms = 15, theta = c(-2, .1), fup = 3,
                        beta = c(2, -0.1, 0.1, -0.2), gamma = 0.6, zeta= c(0.0, -0.2),
                        D = matrix(c(0.25, 0, 0, 0), 2, 2)),
                   simplify = F)

# Now fit using poisson
setwd('../../Multi-test/')
source('EM.R')
poisson.fits <- vector('list', N)
long.formula <- list(Y ~ time + cont + bin + (1|id))
surv.formula <- Surv(survtime, status) ~ bin

for(n in 1:N){
  d <- tests[[n]]$data
  fit <- tryCatch(suppressMessages(EM(long.formula, surv.formula, d, list('poisson'))),
                  error = function(e) NULL)
  poisson.fits[[n]] <- fit
  cat(sprintf('%d/%d\r', n, N))
}

#' Plots --> --------------------------------------------------------------
qz <- qnorm(.975)
# gammas -->
gammas <- sapply(poisson.fits, function(x){
  g <- unname(x$coeff$gamma)
  s <- unname(x$SE['gamma_1'])
  c('lower' = g - qz * s, 'gamma' = g, 'upper' = g + qz * s)
})

plot(gammas[2,], xaxt='n', xlab = '', ylab = bquote(gamma),
     pch = 20, main = bquote(gamma == .(.6)),
     ylim = c(min(gammas), max(gammas)))
abline(h = 0.6, col = 'red', lty = 5)
arrows(x0 = 1:N, 
       y0 = gammas[1,], y1 = gammas[3,],
       angle = 90, code = 3, length = 5e-2)

#' betas -->
betas <- sapply(poisson.fits, function(x){
  b <- unname(x$coeff$beta)
  b
})
SEb <- sapply(poisson.fits, function(x){
  s <- x$SE
  unname(s[grepl('^Y\\_', names(s))])
})

par(mfrow=c(2, 2))
  true.beta <- c(2, -0.1, 0.1, -0.2)
  for(i in 1:4){
    ub <- betas[i,] + qz * SEb[i,]
    lb <- betas[i,] - qz * SEb[i,]
    ii <- i-1
    plot(betas[i,], xaxt='n', xlab = '', ylab = bquote(beta[.(ii)]),
         pch = 20, 
         ylim = c(min(lb), max(ub)))
    abline(h = true.beta[i], col = 'red', lty = 5)
    arrows(x0 = 1:N, 
           y0 = lb, y1 = ub,
           angle = 90, code = 3, length = 5e-2)
  }
par(mfrow=c(1, 1))
