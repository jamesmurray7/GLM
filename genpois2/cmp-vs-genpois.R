rm(list=ls())
assign('gp', get(load('/data/c0061461/GP.RData')))
assign('cmp', get(load('/data/c0061461/CMP.RData')))

lf <- Y ~ time + cont + bin + (1|id)
fit.gp <- glmmTMB(lf, gp, genpois)
fit.cmp <- glmmTMB(lf, cmp, genpois)

# gp is fit with true value -0.5; cmp is fit with true value .5
# gp first
(phi <- exp(fixef(fit.gp)$disp/2)-1)
(delta <- -fixef(fit.cmp)$disp)

-2 * log(1 + phi) # Same level of dispersion in CMP at this level of phi
exp(delta/2)-1

d <- simData_joint2(n = 250, ntms = 10, fup = 3, beta = c(2, -0.1, 0.1, -0.2), delta = c(1.402469 , 0), D = matrix(c(0.25, 0, 0, 0), 2, 2), gamma = 0.6, 
               zeta = c(0, -0.2), theta = c(-2, 0.1), cens.rate = exp(-3.5))$data

fit.cmp.new <-  glmmTMB(lf, d, genpois)
exp(fixef(fit.cmp.new)$disp/2)-1 # appx. -.5, as simulated!


# Do we fit this and get phi = -0.5?---------------------------------------

setwd('../../genpois/')
source('EM.R')
fit.gp <- EM(lf, surv.formula = Surv(survtime, status) ~ bin, data = d,
             control = list(verbose = T, phi.update.quad = T, beta.update.quad = F, include.all = T),
             genpois.inits = T)
a <- my.summary
setwd('../mpcmp/globals/')
source('EM.R')
b <- my.summary
fit.cmp <- EM(lf, disp.formula = ~1, surv.formula = Surv(survtime, status) ~ bin,
              data = d, genpois.inits = T, control = list(verbose = T),
              optim.control = list(optimiser = 'optim', Hessian = 'grad', eps = 1e-3))

sink('/data/c0061461/GP-cmp-compare.txt')
a(fit.gp)
b(fit.cmp)
sink()


# Vice versa - can we simulate on phi = -0.5 and get correct  -------------
# CMP?
setwd('../../genpois/')
source('EM.R')
d <- simData_joint(phi = -0.5)$data
fit.gp <- EM(lf, surv.formula = Surv(survtime, status) ~ bin, data = d,
             control = list(verbose = T, phi.update.quad = T, beta.update.quad = F, include.all = T),
             genpois.inits = T)
setwd('../mpcmp/globals/')
source('EM.R')
b <- my.summary
fit.cmp <- EM(lf, disp.formula = ~1, surv.formula = Surv(survtime, status) ~ bin,
              data = d, genpois.inits = T, control = list(verbose = T),
              optim.control = list(optimiser = 'optim', Hessian = 'grad', eps = 1e-3))
sink('/data/c0061461/GP-cmp-compare2.txt')
a(fit.gp)
b(fit.cmp)
sink()


#' ######################################
#' Full sim study                       #
#' ######################################
rm(list=ls())
save.dir <- '/data/c0061461/18oct-gp-vs-cmp/'
setwd('~/Documents/GLMM/genpois')
source('EM.R')
data <- replicate(100, simData_joint(phi = -0.5)$data, simplify = F)
save(data,  file = paste0(save.dir, 'data.RData'))

long.formula <- Y~time+cont+bin+(1|id)
surv.formula <- Surv(survtime, status) ~ bin
disp.formula <- ~1
genpois.inits <- T
control <- list(verbose = F, phi.update.quad = T, include.all = T)

#' GP ---------------------------------------------------------------
fit.GP <- fit.CMP <- vector('list', 100)
pb <- utils::txtProgressBar(max = 100, style = 3)
for(i in 1:100){
  d <- data[[i]]
  ff <- tryCatch(suppressMessages(EM(long.formula , disp.formula, surv.formula, 
                                     d, control = control, genpois.inits = T)),
                 error = function(e) NULL)
  fit.GP[[i]] <- ff
  utils::setTxtProgressBar(pb, i)
}
close(pb)
save(fit.GP, file = paste0(save.dir, 'GP.RData'))
cat('\n')
#' MPCMP ------------------------------------------------------------
setwd('../mpcmp/globals/')
source('EM.R')
pb <- utils::txtProgressBar(max = 100, style = 3)
for(i in 1:100){
  d <- data[[i]]
  ff <- tryCatch(suppressMessages(EM(long.formula , disp.formula, surv.formula, 
                                     d, control = list(verbose = F, beta.update.quad=F), 
                                     genpois.inits = T,
                                     optim.control = list(optimiser = 'optim', Hessian = 'grad', eps = 1e-3),
                                     summax.fn = function(y) 2 * max(y))),
                 error = function(e) NULL)
  fit.CMP[[i]] <- ff
  utils::setTxtProgressBar(pb, i)
}
close(pb)
save(fit.CMP, file = paste0(save.dir, 'CMP.RData'))
cat('\n')





