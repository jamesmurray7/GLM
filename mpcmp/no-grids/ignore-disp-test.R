rm(list=ls())
source('EM.R')

# Global formulae
disp.formula <- ~1
surv.formula <- Surv(survtime, status) ~ bin

# A function to take a set of data and argument for delta.method (uniroot/optim) and return it

tempfn <- function(data, long.formula, minlength){
  formulas <- parseFormula(long.formula)
  surv <- parseCoxph(surv.formula, data)
  n <- surv$n
  
  #' Initial conditions ----
  inits.long <- Longit.inits(long.formula, disp.formula, data)
  inits.surv <- TimeVarCox(data, inits.long$b, surv$ph, formulas)
  
  # Longitudinal parameters
  beta <- inits.long$beta.init
  D <- inits.long$D.init
  b <- lapply(1:n, function(i) inits.long$b[i, ])
  
  # Survival parameters
  zeta <- inits.surv$inits[match(colnames(surv$ph$x), names(inits.surv$inits))]
  names(zeta) <- paste0('zeta_', names(zeta))
  gamma <- inits.surv$inits[grepl('gamma', names(inits.surv$inits))]
  
  #' Data objects ----
  sv <- surv.mod(surv$ph, surv$survdata, formulas, inits.surv$l0.init)
  dmats <- createDataMatrices(data, formulas, disp.formula)
  # Truncation amount.
  summax <- max(sapply(dmats$Y, max)) * 2
  
  # message('Doing uniroot...')
  # uniroot.inits <- suppressMessages(get.delta.inits(dmats, beta, b, 'uniroot', summax, verbose = T, min.profile.length = 1))
  
  message('Doing optim...')
  optim.inits <- suppressMessages(get.delta.inits(dmats, beta, b, 'optim', summax, verbose = T, min.profile.length = minlength))
  
  # list(
  #   # uni = uniroot.inits,
  #   opt = optim.inits
  # )
  optim.inits
}

# I: Larger beta coeffs --------------------------------------------------
# Heavy under-dispersion
heavy.under <- replicate(100, simData_joint2(n = 250, delta = c(0.8, 0), 
                                            ntms = 10, theta = c(-2, .1), fup = 3,
                                            beta = c(2.0, -0.1, 0.1, -0.2), gamma = 0.6, zeta= c(0.0, -0.2),
                                            D = matrix(c(0.25, 0, 0, 0.00), 2, 2))$data,
                         simplify = F)

heavy.under.estimates <- lapply(heavy.under, tempfn, Y ~ time + cont + bin + (1|id), 1)

# light under-dispersion
light.under <- replicate(100, simData_joint2(n = 250, delta = c(0.3, 0), 
                                            ntms = 10, theta = c(-2, .1), fup = 3,
                                            beta = c(2.0, -0.1, 0.1, -0.2), gamma = 0.6, zeta= c(0.0, -0.2),
                                            D = matrix(c(0.25, 0, 0, 0.00), 2, 2))$data,
                         simplify = F)

light.under.estimates <- lapply(light.under, tempfn, Y ~ time + cont + bin + (1|id), 1)

# light over-dispersion
light.over <- replicate(100, simData_joint2(n = 250, delta = c(-0.2, 0), 
                                            ntms = 10, theta = c(-2, .1), fup = 3,
                                            beta = c(2.0, -0.1, 0.1, -0.2), gamma = 0.6, zeta= c(0.0, -0.2),
                                            D = matrix(c(0.25, 0, 0, 0.00), 2, 2))$data,
                         simplify = F)

light.over.estimates <- lapply(light.over, tempfn, Y ~ time + cont + bin + (1|id), 1)

# Report median/mean with/without the inclusion of UB/LB estimates...
report <- function(x){
  cat(sprintf('The median value for all estimates is %.3f and mean value is %.3f.\n', x$media, x$mean))
  o <- x$subjec
  to.keep <- abs(round(o, 3)) < 2 & !is.na(o)
  o <- o[to.keep]
  cat(sprintf('Using only the %d/%d estimates not at the lower/upper bounds (+/- 2)\nWe obtain median %.3f and mean %.3f\n----\n',
              length(o), length(x$subject), median(o), mean(o)))
}

invisible(lapply(heavy.under.estimates, report)) # True value: 0.8
invisible(lapply(light.under.estimates, report)) # True value: 0.3
invisible(lapply(light.over.estimates, report))  # True value: -0.2

# See how close estimates are (return some matrix)
how.close <- function(x, target){
  # median and mean of ALL returned data
  o <- x$subjec
  to.keep <- abs(round(o, 3)) < 2 & !is.na(o)
  o <- o[to.keep]
  data.frame(median_all = x$median.estimate, 
             median_cut = median(o, na.rm = T),
             mean_all = x$mean.estimate,
             mean_cut = mean(o, na.rm = T)) - target
}

aa <- do.call(rbind, lapply(heavy.under.estimates, how.close, .8))
bb <- do.call(rbind, lapply(light.under.estimates, how.close, .3))
cc <- do.call(rbind, lapply(light.over.estimates, how.close, -.2))

colMeans(aa); colMeans(bb); colMeans(cc) # Looks like median_cut generates best swathe of estimates.

#' =============================
# Heavier overdispersion...
heavy.over <- replicate(10, simData_joint2(n = 250, delta = c(-0.5, 0), 
                                           ntms = 10, theta = c(-2, .1), fup = 3,
                                           beta = c(3.0, -0.1, 0.1, -0.2), gamma = 0.6, zeta= c(0.0, -0.2),
                                           D = matrix(c(0.25, 0, 0, 0.00), 2, 2))$data,
                        simplify = F)

heavy.over.estimates <- lapply(heavy.over, tempfn, Y ~ time + cont + bin + (1|id), 1)
invisible(lapply(heavy.over.estimates, report))
do.call(rbind, lapply(heavy.over.estimates, how.close, -0.5))

#' =======================
#' Works quite well with _any_ large(r) values, slowed down largely by summax amount.
#' =======================

# 24/8/22 functions to plot -----------------------------------------------
plot.density.ests <- function(l, cut = F, true.delta){ # l is a list.
  ests <- lapply(l, el)
  if(cut) 
    ests <- lapply(ests, function(x) x[round(abs(x), 3) < 2 & !is.na(x)])
  
  # Work out densities of subject-specific point estimates...
  densities <- lapply(ests, density)
  
  # Limits 
  xlims <- c(min(sapply(densities, function(d) min(d$x))),
             max(sapply(densities, function(d) max(d$x))))
  ylims <- c(min(sapply(densities, function(d) min(d$y))),
             max(sapply(densities, function(d) max(d$y))))
  
  # Plot
  title <- bquote("True " ~ delta == .(true.delta))
  plot(density(ests[[1]]), xlab = expression(delta), main = title, 
       xlim = xlims, ylim = ylims)
  for(i in 2:length(densities))
    lines(densities[[i]], col = i)
  abline(v = true.delta)
}

# hu: heavy under; lu: light under; lo: light over.
times.hu <- sapply(heavy.under.estimates, function(x) x$time)
times.lu <- sapply(light.under.estimates, function(x) x$time)
times.lo <- sapply(light.over.estimates, function(x) x$time)

boxplot(times.hu, times.lu, times.lo,
        ylab = 'seconds', main = 'Elapsed time', xaxt = 'n')
axis(1, at = c(1,2,3), labels = c(expression(delta==0.8),expression(delta==0.3),expression(delta==-0.2)))

