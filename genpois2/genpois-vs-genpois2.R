rm(list=ls())
assign('gp1', get(load('gp1set.RData'))) # BE IN GENPOIS(1)
assign('tmb', get(load('../genpois2/gp1set.RData')))
assign('cmp', get(load('../mpcmp/globals/mpcmpset.RData')))

2*log(gp1[[1]]$coeffs$phi+1)       # phi -> delta;
exp(tmb[[1]]$coeffs$delta/2)-1     # delta -> phi.




# Plots -------------------------------------------------------------------
targets <- data.frame(target = c(.25, 2, -0.1, 0.1, -0.2, -0.5, 0.6, -0.2),
                      parameter = names(tmb[[1]]$SE))

ests.tmb <- as.data.frame(do.call(rbind, lapply(tmb, function(x){
  co <- x$coeff
  setNames(c(co$D, c(co$beta), exp(co$delta/2)-1, co$gamma, co$zeta),
           names(x$SE))
})))
ests.tmb$method <- 'TMB'

ests.gp1 <- as.data.frame(do.call(rbind, lapply(gp1, function(x){
  co <- x$coeff
  setNames(c(co$D, c(co$beta), co$phi, co$gamma, co$zeta),
           c(names(x$SE)[1:5], 'delta_(Intercept)', names(x$SE)[7:8]))
})))
ests.gp1$method <- 'GP-1'

ests.cmp <- as.data.frame(do.call(rbind, lapply(cmp, function(x){
  co <- x$coeff
  setNames(c(co$D, c(co$beta), exp(-co$delta/2)-1, co$gamma, co$zeta),
           c(names(x$SE)[1:5], 'delta_(Intercept)', names(x$SE)[7:8]))
})))
ests.cmp$method <- 'MPCMP'

library(tidyverse)
all.ests <- rbind(ests.tmb, ests.gp1, ests.cmp)
all.ests %>% 
  pivot_longer(`D[1,1]`:`zeta_bin`, names_to = 'parameter', values_to = 'estimate') %>% 
  left_join(., targets, 'parameter') %>% 
  ggplot(aes(x = estimate, colour = method)) + 
  geom_vline(aes(xintercept = target)) + 
  geom_density() + 
  facet_wrap(~parameter, scales= 'free')+
  theme_light()
  
# Timings -----------------------------------------------------------------
et.tmb <- as.data.frame(do.call(rbind, lapply(tmb, function(x) x$elapsed.time)))
et.tmb$method <- 'TMB'
et.gp1 <- as.data.frame(do.call(rbind, lapply(gp1, function(x) x$elapsed.time)))
et.gp1$method <- 'GP-1'

timings <- rbind(et.tmb, et.gp1)
timings %>% 
  pivot_longer(`startup time`:`Total computation time`, 
               names_to = 'Timing', values_to = 'elapsed') %>% 
  filter(Timing %in% c('EM time', 'Total computation time')) %>% 
  ggplot(aes(x = Timing, fill = method, y = elapsed)) + 
  geom_boxplot()


# Coverage ----------------------------------------------------------------
qz <- qnorm(.975)
# TMB
colSums(
  do.call(rbind, lapply(tmb, function(x){
    co <- x$coeff; se <- x$SE
    all.but.phi <- c(co$D, c(co$beta), co$gamma, co$zeta)
    all.but.phi.lb <- all.but.phi - qz * x$SE[!grepl('^delta', names(x$SE))]
    all.but.phi.ub <- all.but.phi + qz * x$SE[!grepl('^delta', names(x$SE))]
    # phis: 
    del <- co$delta
    lb <- del - qz * x$SE[grepl('^delta', names(x$SE))]
    ub <- del + qz * x$SE[grepl('^delta', names(x$SE))]
    
    lb <- exp(lb/2)-1; ub <- exp(ub/2)-1
    
    lbs <- setNames(c(all.but.phi.lb[1:5], lb, all.but.phi.lb[6:7]),
                    names(se))
    ubs <- setNames(c(all.but.phi.ub[1:5], ub, all.but.phi.ub[6:7]),
                    names(se))
    lbs <= targets$target & ubs >= targets$target
  }))
)/100

# GP-1
colSums(
  do.call(rbind, lapply(gp1, function(x){
    co <- x$coeff; se <- x$SE
    est <- setNames(c(co$D, c(co$beta), co$phi, co$gamma, co$zeta),
             c(names(x$SE)[1:5], 'delta_(Intercept)', names(x$SE)[7:8]))
    
    lb <- est - qz * se; ub <- est + qz * se
    
    lb <= targets$target & ub >= targets$target
  }))
)/100

