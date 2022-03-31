rm(list=ls())
source('EM.R')
diag(true.D) <- c(.50^2, .05, .3^2, .05, .50^2, .05)
true.D <- as.matrix(Matrix::nearPD(true.D)$mat)

# Load the two fits(N)Q RData files
load('~/Downloads/Quad-fits-29mar.RData')
load('~/Downloads/Noquad-fits-29mar.RData')

extract.coeffs <- function(x){ # x a sub-list of one model fit
  co <- x$co
  setNames(c(vech(co$D), c(co$beta), co$var.e, co$gamma, co$eta),
           names(x$SE))
}
extract.SE <- function(x) x$SE # x a sub-list of one model fit

Q.ests <- lapply(fitsQ, function(f) do.call(rbind, lapply(f, extract.coeffs)))
Q.SE <- lapply(fitsQ, function(f) do.call(rbind, lapply(f, extract.SE)))                 
NQ.ests <- lapply(fitsNQ, function(f) do.call(rbind, lapply(f, extract.coeffs)))
NQ.SE <- lapply(fitsNQ, function(f) do.call(rbind, lapply(f, extract.SE)))                 

targets <- apply(matrix(c(vech(true.D), do.call(c, lapply(1:3, function(i) true.beta[i,])),
                          0.25, true.gamma, true.eta),nr=1),
                 2,rep,100)

tabulate <- function(tab, SE, target.col = T){
  m <- apply(tab,2,mean)
  SD <- apply(tab,2,sd)
  SEs <- apply(SE,2,mean)
  ub <- tab + qnorm(.975) * SE; lb <- tab - qnorm(.975) * SE
  CP <- colSums(ub >= targets & lb <= targets)/100
  out <- cbind(mean = m, SD = SD, SE = SEs, CP = CP)
  if(target.col) out <- cbind(target = targets[1,], out)
  out
}

quad_20pc <- tabulate(Q.ests[[1]], Q.SE[[1]], T)
quad_40pc <- tabulate(Q.ests[[2]], Q.SE[[2]], T)
noquad_20pc <- tabulate(NQ.ests[[1]], NQ.SE[[1]], F)
noquad_40pc <- tabulate(NQ.ests[[2]], NQ.SE[[2]], F)

pc20 <- cbind(quad_20pc, noquad_20pc)
pc40 <- cbind(quad_40pc, noquad_40pc)

library(tidyverse)

colnames(pc20)[6:9] <- paste0('n', colnames(pc20)[6:9])
colnames(pc40)[6:9] <- paste0('n', colnames(pc40)[6:9])

todp <- function(x, n) format(round(x, n), nsmall = n)

pc20df <- cbind(parameter = rownames(pc20), pc20 %>% 
  as.data.frame %>% 
  mutate(across(-contains('CP'), ~ todp(.x, 3)),
         across(contains('CP'),  ~ todp(.x, 2)))) %>% 
  as_tibble %>%
  mutate(
    parameter = gsub('^cont$', 'zeta_1', parameter),
    parameter = gsub('^bin$', 'zeta_2', parameter),
    parameter = gsub('\\(Intercept\\)', '0}', parameter),
    parameter = gsub('time$', '1}', parameter),
    parameter = gsub('cont$', '2}', parameter),
    parameter = gsub('bin$', '3}', parameter),
    parameter = gsub('\\[', '_{', parameter),
    parameter = gsub('\\]', '}', parameter),
    parameter = gsub('^G\\_', 'beta_{1', parameter),
    parameter = gsub('^B\\_', 'beta_{2', parameter),
    parameter = gsub('^P\\_', 'beta_{3', parameter),
    parameter = gsub('var\\.e', 'sigma^2_varepsilon', parameter),
    parameter = paste0('$', parameter, '$')
  )

pc40df <- cbind(parameter = rownames(pc40), pc40 %>% 
                  as.data.frame %>% 
                  mutate(across(-contains('CP'), ~ todp(.x, 3)),
                         across(contains('CP'),  ~ todp(.x, 2)))) %>% 
  as_tibble %>%
  mutate(
    parameter = gsub('^cont$', 'zeta_1', parameter),
    parameter = gsub('^bin$', 'zeta_2', parameter),
    parameter = gsub('\\(Intercept\\)', '0}', parameter),
    parameter = gsub('time$', '1}', parameter),
    parameter = gsub('cont$', '2}', parameter),
    parameter = gsub('bin$', '3}', parameter),
    parameter = gsub('\\[', '_{', parameter),
    parameter = gsub('\\]', '}', parameter),
    parameter = gsub('^G\\_', 'beta_{1', parameter),
    parameter = gsub('^B\\_', 'beta_{2', parameter),
    parameter = gsub('^P\\_', 'beta_{3', parameter),
    parameter = gsub('var\\.e', 'sigma^2_varepsilon', parameter),
    parameter = paste0('$', parameter, '$')
  )

library(xtable)
x20 <-xtable(pc20df); x40 <- xtable(pc40df)
print(x20, 
      include.rownames = F,
      sanitize.text.function = identity)

print(x40, 
      include.rownames = F,
      sanitize.text.function = identity)
