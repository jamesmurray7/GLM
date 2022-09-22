load('~/Downloads/intonly-poiss-vs-mpcmp-poissfit.RData')
load('~/Downloads/intonly-poiss-vs-mpcmp-mpcmpfit.RData')
load('~/Downloads/intonly-poiss-vs-mpcmp-mpcmpfit2.RData') # This one has unique summaxes
                                                             # (i.e. ../half_and_half2/*)

#' Quick overviews --------------------------------------------------------
# How many are missing?
sum(sapply(poisson.fit, is.null)) / 50  # All seem relatively stable methods, then.
sum(sapply(mpcmp.fit, is.null))   / 50      
sum(sapply(mpcmp2.fit, is.null))  / 50

#' Average Estimate and SE for \gamma -------------------------------------
g.poisson <- sapply(poisson.fit, function(x) unname(x$coeff$gamma))
g.mpcmp <- sapply(mpcmp.fit, function(x) unname(x$coeff$gamma))
g.mpcmp2 <- sapply(mpcmp2.fit, function(x) unname(x$coeff$gamma))
plot(g.poisson, g.mpcmp, pch = 20, main = bquote(hat(gamma)),
     xlab = 'Poisson', ylab = 'MPCMP'); abline(0, 1) # MPCMP nearly always bigger.
plot(g.mpcmp, g.mpcmp2 , pch = 20, main = bquote(hat(gamma)),
     xlab = 'MPCMP', ylab = 'MPCMP2'); abline(0, 1) # MPCMP nearly always bigger.
# Standard error
g.se.poisson <- sapply(poisson.fit, function(x) unname(x$SE['gamma_1']))
g.se.mpcmp <- sapply(mpcmp.fit, function(x) unname(x$SE['gamma']))
g.se.mpcmp2 <- sapply(mpcmp2.fit, function(x) unname(x$SE['gamma']))
plot(g.se.poisson, g.se.mpcmp, pch = 20, main = bquote('SE['*gamma*']'),
     xlab = 'Poisson', ylab = 'MPCMP'); abline(0, 1) # MPCMP nearly always bigger.
plot(g.se.mpcmp,  g.se.mpcmp2, pch = 20, main = bquote('SE['*gamma*']'),
     xlab = 'MPCMP', ylab = 'MPCMP2'); abline(0, 1) 

# Overall (95%) coverage
qz <- qnorm(.975)
target <- c(.25, 2., -.1, .1, -.2, .6, -.2)
mpcmp <- lapply(mpcmp.fit, function(x){
  co <- x$coeff
  se <- x$SE
  ests <- setNames(
    c(co$D, c(co$beta), co$gamma, co$zeta),
    names(se)
  )
  lb <- ests - qz * se; ub <- ests + qz * se
  lb <= target & ub >= target
})

mpcmp2 <- lapply(mpcmp2.fit, function(x){
  co <- x$coeff
  se <- x$SE
  ests <- setNames(
    c(co$D, c(co$beta), co$gamma, co$zeta),
    names(se)
  )
  lb <- ests - qz * se; ub <- ests + qz * se
  lb <= target & ub >= target
})

poiss <- lapply(poisson.fit, function(x){
  co <- x$coeff
  se <- x$SE
  ests <- setNames(
    c(co$D, c(co$beta), co$gamma, co$zeta),
    names(se)
  )
  lb <- ests - qz * se; ub <- ests + qz * se
  lb <= target & ub >= target
})

apply(do.call(rbind, poiss),  2, sum)/50 * 100
apply(do.call(rbind, mpcmp),  2, sum)/50 * 100
apply(do.call(rbind, mpcmp2), 2, sum)/50 * 100

#' Plots --> --------------------------------------------------------------
#' \gammas -->
p.gammas <- sapply(poisson.fit, function(x){
  g <- unname(x$coeff$gamma)
  s <- unname(x$SE['gamma_1'])
  c('lower' = g - qz * s, 'gamma' = g, 'upper' = g + qz * s)
})

cmp.gammas <- sapply(mpcmp.fit, function(x){
  g <- unname(x$coeff$gamma)
  s <- unname(x$SE['gamma'])
  c('lower' = g - qz * s, 'gamma' = g, 'upper' = g + qz * s)
})

cmp2.gammas <- sapply(mpcmp2.fit, function(x){
  g <- unname(x$coeff$gamma)
  s <- unname(x$SE['gamma'])
  c('lower' = g - qz * s, 'gamma' = g, 'upper' = g + qz * s)
})

par(mfrow=c(2,1), mai = c(.45, 1, .5, 1))
  # poisson gammas
  plot(p.gammas[2,], xaxt='n', xlab = '', ylab = bquote(gamma),
       pch = 20,
       ylim = c(min(p.gammas), max(p.gammas)))
  title(main = bquote(gamma==.(.6)), xlab = 'Poisson', mgp = c(1,1,0))
  abline(h = 0.6, col = 'red', lty = 5)
  arrows(x0 = 1:50, 
         y0 = p.gammas[1,], y1 = p.gammas[3,],
         angle = 90, code = 3, length = 5e-2)
  # MPCMP(1) gammas
  plot(cmp.gammas[2,], xaxt='n', xlab = '', ylab = bquote(gamma),
       pch = 20, main = '',
       ylim = c(min(cmp.gammas), max(cmp.gammas)))
  title(xlab = 'MPCMP', mgp = c(1,1,0))
  abline(h = 0.6, col = 'red', lty = 5)
  arrows(x0 = 1:50, 
         y0 = cmp.gammas[1,], y1 = cmp.gammas[3,],
         angle = 90, code = 3, length = 5e-2)
par(mfrow=c(1,1))

#' \betas -->
# Poisson
p.betas <- sapply(poisson.fit, function(x){
  b <- unname(x$coeff$beta)
  b
})
p.SEb <- sapply(poisson.fit, function(x){
  s <- x$SE
  unname(s[grepl('^Y\\_', names(s))])
})
# MPCMP
m.betas <- sapply(mpcmp.fit, function(x){
  b <- unname(x$coeff$beta)
  b
})
m.SEb <- sapply(mpcmp.fit, function(x){
  s <- x$SE
  unname(s[grepl('^beta\\_', names(s))])
})

# Trying to plot all...
par(mfrow=c(2, 2))
true.beta <- c(2, -0.1, 0.1, -0.2)
for(i in 1:4){
  p.ub <- p.betas[i,] + qz * p.SEb[i,]; m.ub <- m.betas[i,] + qz * m.SEb[i,]
  p.lb <- p.betas[i,] - qz * p.SEb[i,]; m.lb <- m.betas[i,] - qz * m.SEb[i,]
  ii <- i-1
  plot(p.betas[i,], xaxt='n', xlab = '', ylab = bquote(beta[.(ii)]),
       pch = 20, 
       ylim = c(min(pmin(p.lb, m.lb)), max(pmax(p.ub, m.ub))))
  points(1:50 + 5e-1, m.betas[i,], pch = 20, col = 'blue')
  abline(h = true.beta[i], col = 'red', lty = 5)
  arrows(x0 = 1:50, 
         y0 = p.lb, y1 = p.ub,
         angle = 90, code = 3, length = 5e-2)
  arrows(x0 = 1:50 + 5e-1, 
         y0 = m.lb, y1 = m.ub,
         angle = 90, code = 3, length = 5e-2, col = 'blue')
} # A bit messy!
par(mfrow=c(1, 1))
