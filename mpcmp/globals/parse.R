#' ####
#' Parse-sim-study.R
#' Parsing results from sim-study.R
#' ####

rm(list=ls())
data.dir <- '/data/c0061461/14oct22/'
load(paste0(data.dir, 'global_intslope.RData'))


# Plot parameter estimates ------------------------------------------------
library(dplyr)
library(ggplot2)
library(tidyr)

# Set out true values
true.value <- data.frame(
  parameter = names(fits[[1]]$SE),
  true = c(.25, 2, -0.1, 0.1, -0.2, .6, -.1, .6, -.2)
)

# Collate all fitted objects...
ests <- as.data.frame(do.call(rbind, lapply(fits, function(x){
  co <- x$coeffs
  setNames(c(co$D, c(co$beta), c(co$delta), co$gamma, co$zeta),
           names(x$SE))
})))

# Joining/cleaning.
ests <- ests %>% 
  pivot_longer(`D[1,1]`:`zeta_bin`, names_to = 'parameter', values_to = 'estimate') %>% 
  left_join(., true.value, 'parameter')

#' The plot
ests %>% 
  ggplot(aes(x = estimate)) + 
  geom_vline(aes(xintercept = true)) + 
  geom_density() + 
  facet_wrap(~parameter, scales = 'free') + 
  theme_bw()

ggsave('/data/c0061461/sept-29-sims/parameter_estimates.png')

# Parameter coverage ------------------------------------------------------
qz <- qnorm(.975)
CP <- colSums(do.call(rbind, lapply(fits, function(x){
  co <- x$coeffs; se <- x$SE
  est <- setNames(c(co$D, c(co$beta), c(co$delta), co$gamma, co$zeta),
           names(se))
  lb <- est - qz * se; ub <- est + qz * se
  lb <= true.value$true & ub >= true.value$true
})))/100

CP # 94% coverage on deltas; why is beta_intercept flagging?

