rm(list=ls())
source('EM.R')

# Load the two fits(N)Q RData files
load('~/Downloads/fits-rustand-failure2040.RData')


extract.coeffs <- function(x){ # x a sub-list of one model fit
  co <- x$co
  setNames(c(vech(co$D), c(co$beta), co$var.e, co$gamma),
           names(x$SE))
}
extract.SE <- function(x) x$SE # x a sub-list of one model fit

pc40.ests <- do.call(rbind, lapply(fits$`40%`, extract.coeffs))
pc40.SE <-   do.call(rbind, lapply(fits$`40%`, extract.SE))            
pc20.ests <- do.call(rbind, lapply(fits$`20%`, extract.coeffs))
pc20.SE <-   do.call(rbind, lapply(fits$`20%`, extract.SE))                         

targets <- apply(matrix(c(vech(true.D), do.call(c, lapply(1:3, function(i) true.beta[i,])),
                          0.16, true.gamma),nr=1),
                 2,rep,100)

tabulate <- function(tab, SE, target.col = T){
  m <- apply(tab,2,mean)
  SD <- apply(tab,2,sd)
  SEs <- apply(SE,2,mean)
  ub <- tab + qnorm(.975) * SE; lb <- tab - qnorm(.975) * SE
  CP <- colSums(ub >= targets & lb <= targets)/100
  # biases
  bias <- tab - targets
  mean.bias <- apply(bias, 2, mean)
  sd.bias <- apply(bias, 2, sd)
  # form output
  out <- cbind(mean = m, SD = SD, SE = SEs, mbias = mean.bias, sbias = sd.bias, CP = CP)
  if(target.col) out <- cbind(target = targets[1,], out)
  out
}

pc40 <- tabulate(pc40.ests, pc40.SE, T)
pc20 <- tabulate(pc20.ests, pc20.SE, T)

library(tidyverse)
todp <- function(x, n) format(round(x, n), nsmall = n)

tidy.out <- function(tab){
  if(!"data.frame"%in%class(tab)) x <- as.data.frame(tab)
  if(is.null(as.data.frame(tab)$parameter)) x$parameter <- rownames(tab)
  x
  x$paramtarget <- paste0(x$parameter, ' = ', x$target)
  x <- x %>% 
    mutate(
      across(mean:sbias, ~ todp(.x, 3)),
      CP = todp(CP, 2),
      paramtarget = paste0(parameter, ' = ', target),
      bias = paste0(mbias, ' (', sbias, ')')
    ) %>% 
    select(-parameter, -target) %>% 
    select(parameter = paramtarget, Mean = mean, SD = SD, SE = SE, Bias = bias, CP = CP)
  
  # Fixing parameter
  x <- x %>% 
    mutate(
      parameter = gsub('\\(Intercept\\)', '0}', parameter),
      parameter = gsub('time', '1}', parameter),
      parameter = gsub('cont', '2}', parameter),
      parameter = gsub('bin', '3}', parameter),
      parameter = gsub('\\[', '_{', parameter),
      parameter = gsub('\\]', '}', parameter),
      parameter = gsub('^G\\_', 'beta_{1', parameter),
      parameter = gsub('^B\\_', 'beta_{2', parameter),
      parameter = gsub('^P\\_', 'beta_{3', parameter),
      parameter = gsub('var\\.e', 'sigma^2_varepsilon', parameter),
      parameter = paste0('$', parameter, '$')
    )
  
  return(x)
}

pc20df <- tidy.out(pc20)
pc40df <- tidy.out(pc40)

df <- left_join(pc20df, pc40df, 'parameter')

library(xtable)

print(xtable(df), 
      include.rownames = F,
      sanitize.text.function = identity)

print(x40, 
      include.rownames = F,
      sanitize.text.function = identity)
