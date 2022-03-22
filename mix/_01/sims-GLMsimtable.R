# Tabulate Sims
rm(list=ls())
load('~/Downloads/mixfits-temp-40pc.RData')
source('EM.R')
diag(true.D) <- c(.50^2, .05, .3^2, .05, .50^2, .05)
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

names(fits) <- c('n=250', 'n=500')#, 'n=1000')

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
for(i in 1:2){
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

library(tidyverse)
dp <- function(x, n) format(round(x, n), nsmall = n)

fit.table2 <- lapply(fit.table, function(tab){
  id <- gsub('n\\=', '', tab$id)
  tab <- tab %>% mutate(across(mean:CI.up, ~ dp(.x, 3)))
  tab <- tab %>% mutate(CP =dp(CP, 2))
  tab$CI = paste0('[', tab$CI.low, ', ', tab$CI.up, ']')
  tab %>% 
    select(parameter, Mean = mean, SD = SD, SE = SE, CI = CI, CP = CP)
})

fit.table3 <- lapply(fit.table2, function(tab){
  tab$parameter <- gsub('\\[', '_{', tab$parameter)
  tab$parameter <- gsub('\\]', '}', tab$parameter)
  tab$parameter <- gsub('^cont$', 'zeta_1', tab$parameter)
  tab$parameter <- gsub('^bin$', 'zeta_2', tab$parameter)
  tab$parameter <- gsub('^G\\_', 'beta_{1', tab$parameter)
  tab$parameter <- gsub('^B\\_', 'beta_{2', tab$parameter)
  tab$parameter <- gsub('^P\\_', 'beta_{3', tab$parameter)
  tab$parameter <- gsub('\\(Intercept\\)', '0}', tab$parameter)
  tab$parameter <- gsub('time', '1}', tab$parameter)
  tab$parameter <- gsub('cont', '2}', tab$parameter)
  tab$parameter <- gsub('bin',  '3}', tab$parameter)
  tab$parameter <- gsub('var\\.e',  'sigma^2_varepsilon', tab$parameter)
  tab$parameter <- gsub('^gamma', 'gamma', tab$parameter)
  tab$parameter <- paste0('$', tab$parameter, '$')
  tab
})

tab1 <- fit.table3[[1]]
tab1 <- cbind(Parameter = tab1[,1], Target = dp(targets,3), tab1[,-1])

tab <- cbind(tab1, fit.table3[[2]][,-1])

library(xtable)
print(xtable(tab),
      include.rownames = F,
      sanitize.text.function = identity)


# Time taken --------------------------------------------------------------

sapply(1:2, function(i){
  f1 <- fits[[i]]
  n1 <- names(fits)[i]
  et <- do.call(c, lapply(f1, function(f) f$EMtime + f$postpro))
  m <- median(et)
  p25 <- unname(quantile(et, .25)) ; p75 <- unname(quantile(et, .75))
  
  cat(sprintf(
    "%s Median [IQR]: %.1f(seconds), [%.1f, %.1f]\n\n", n1, m, p25, p75
  ))
  
})





