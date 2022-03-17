# Tabulate Sims
rm(list=ls())
load('~/Downloads/mixfits-2.RData')
source('EM.R')
diag(true.D) <- c(.25, .06, .50, .04, .25, .05)
true.D <- as.matrix(Matrix::nearPD(true.D)$mat)


# Setting target values (known) -------------------------------------------
pnams <- names(fits[[1]][[1]]$SE)

targets <- setNames(c(
  vech(true.D),
  do.call(c, lapply(1:3, function(i) true.beta[i,])),
  0.25,
  true.gamma, true.eta), pnams)

# 100 possible fits.
targets.mat <- apply(t(as.matrix(targets)), 2, rep, 100)


# Getting fits: Coeffs, SD, SE, CI, CP ------------------------------------

names(fits) <- c('n=250', 'n=500', 'n=1000')

extract.coeffs <- function(x){ # x a 'sub-list'
  xc <- x$coeffs
  vD <- vech(xc$D)
  names(vD) <- paste0('[' , apply(which(lower.tri(xc$D, T), arr.ind = T), 1, paste0, collapse = ', '), ']')
  beta <- xc$beta
  names(beta) <- paste0(rep(c('G_', 'B_', 'P_'), each = 4), c('(Intercept)', 'time', 'cont', 'bin'))
  out <- c(vD, beta, 'var.e' = xc$var.e, xc$gamma, xc$eta)
  out
}

# Can't think how to double lapply, so let's just for loop
fit.table <- list()
for(i in 1:3){
  fit <- fits[[i]]
  ests <- do.call(rbind, 
                  lapply(fit, function(x) if(!is.null(x)) extract.coeffs(x)))
  SE <- do.call(rbind, 
                lapply(fit, function(x) if(!is.null(x)) x$SE))
  
  mean <- apply(ests, 2, mean)
  emp.SD <- apply(ests, 2, sd)
  SEs <- apply(SE, 2, mean)
  lb <- ests - qnorm(.975) * SE; ub <- ests + qnorm(.975) * SE
  CP <- apply(lb <= targets.mat & ub >= targets.mat, 2, sum) / 100
  
  
  cat(sprintf('%d: %d successful fits out of 100\n\n', i, nrow(ests)))
  
  fit.table[[i]] <- data.frame(parameter = colnames(SE), mean = mean, SD = emp.SD, SE = SEs, 
                          CI.low = apply(lb, 2, mean), CI.up = apply(ub, 2, mean),
                          CP = CP, id = names(fits)[i])
}
