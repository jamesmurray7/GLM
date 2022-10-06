#' ####
#' Parse-sim-study.R
#' Parsing results from sim-study.R
#' ####

rm(list=ls())
data.dir <- '/data/c0061461/sept-29-sims/'
dir(data.dir, pattern = '^fit')

#' Check for null fits ---
sum(sapply(fits.inits.quad, is.null))
sum(sapply(fits.none.delta, is.null))
sum(sapply( fits.none.none, is.null))    # All zero, which is good!
sum(sapply(   fits.poisson, is.null))



# Plot parameter estimates ------------------------------------------------
library(dplyr)
library(ggplot2)
library(tidyr)

# Set out true values
true.value <- data.frame(
  parameter = names(fits.inits.quad[[1]]$SE),
  true = c(.25, 2, -0.1, 0.1, -0.2, .6, -.2)
)

# Collate all fitted objects...
ests.inits.quad <- as.data.frame(do.call(rbind, lapply(fits.inits.quad, function(x){
  co <- x$coeffs
  setNames(c(co$D, c(co$beta), co$gamma, co$zeta),
           names(x$SE))
})))

ests.none.delta <- as.data.frame(do.call(rbind, lapply(fits.none.delta, function(x){
  co <- x$coeffs
  setNames(c(co$D, c(co$beta), co$gamma, co$zeta),
           names(x$SE))
})))

ests.none.none <- as.data.frame(do.call(rbind, lapply(fits.none.none, function(x){
  co <- x$coeffs
  setNames(c(co$D, c(co$beta), co$gamma, co$zeta),
           names(x$SE))
})))

ests.poisson <- as.data.frame(do.call(rbind, lapply(fits.poisson, function(x){
  co <- x$coeffs
  setNames(c(co$D, c(co$beta), co$gamma, co$zeta),
           c("D[1,1]", "beta_(Intercept)", "beta_time", "beta_cont", "beta_bin", "gamma", "zeta_bin"))
})))

# Joining/cleaning.
all.ests <- rbind(
  ests.inits.quad %>% mutate(method = 'Inits, double quad'),
  ests.none.delta %>% mutate(method = 'No inits, delta quad'),
  ests.none.none %>% mutate(method = 'No inits, no quad'),
  ests.poisson %>% mutate(method = 'Poisson')
) %>% 
  pivot_longer(`D[1,1]`:`zeta_bin`, names_to = 'parameter', values_to = 'estimate') %>% 
  left_join(., true.value, 'parameter')

#' The plot
all.ests %>% 
  ggplot(aes(x = estimate, colour = method)) + 
  geom_vline(aes(xintercept = true)) + 
  geom_density() + 
  facet_wrap(~parameter, scales = 'free') + 
  theme_bw()

ggsave('/data/c0061461/sept-29-sims/parameter_estimates.png')

# Parameter coverage ------------------------------------------------------
qz <- qnorm(.975)
CP.inits.quad <- colSums(do.call(rbind, lapply(fits.inits.quad, function(x){
  co <- x$coeffs; se <- x$SE
  est <- setNames(c(co$D, c(co$beta), co$gamma, co$zeta),
           names(se))
  lb <- est - qz * se; ub <- est + qz * se
  lb <= true.value$true & ub >= true.value$true
})))/100

CP.none.delta <- colSums(do.call(rbind, lapply(fits.none.delta, function(x){
  co <- x$coeffs; se <- x$SE
  est <- setNames(c(co$D, c(co$beta), co$gamma, co$zeta),
                  names(se))
  lb <- est - qz * se; ub <- est + qz * se
  lb <= true.value$true & ub >= true.value$true
})))/100

CP.none.none <- colSums(do.call(rbind, lapply(fits.none.none, function(x){
  co <- x$coeffs; se <- x$SE
  est <- setNames(c(co$D, c(co$beta), co$gamma, co$zeta),
                  names(se))
  lb <- est - qz * se; ub <- est + qz * se
  lb <= true.value$true & ub >= true.value$true
})))/100

CP.poisson <- colSums(do.call(rbind, lapply(fits.poisson, function(x){
  co <- x$coeffs; se <- x$SE
  est <- setNames(c(co$D, c(co$beta), co$gamma, co$zeta),
                  c("D[1,1]", "beta_(Intercept)", "beta_time", "beta_cont", "beta_bin", "gamma", "zeta_bin"))
  lb <- est - qz * se; ub <- est + qz * se
  lb <= true.value$true & ub >= true.value$true
})))/100

sink('/data/c0061461/sept-29-sims/CPs.txt')
print('CP inits quad')
print(CP.inits.quad)
print('CP.none.delta')
print(CP.none.delta)
print('CP.none.none')
print(CP.none.none)
print('CP.poisson')
print(CP.poisson)
sink()


# Dispersion --------------------------------------------------------------
fitted.disps(fits.inits.quad[[1]], tests[[1]])
all.fits <- list(fits.inits.quad, fits.none.delta, fits.none.none)
 
# Approximately 90% coverage of true dispersion values.
delta.cps <- do.call(cbind, lapply(all.fits, function(L){
  this.list <- lapply(1:100, function(l){
    ff <- invisible(fitted.disps(L[[l]], tests[[l]], T, F))
    sum(ff[ff$truncated == '', 'coverage'])/sum(ff$truncated=='')
  })
  do.call(rbind, this.list)
}))
colnames(delta.cps) <- c('Inits + quad', 'No inits + quad', 'No inits, no quad')

sink('/data/c0061461/sept-29-sims/deltaCP.txt')
delta.cps
sink()

png('/data/c0061461/sept-29-sims/delta_CPs.png')
boxplot(delta.cps,
        xlab = '', ylab = '95% CP')
dev.off()
