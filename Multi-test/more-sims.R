source('simData.R')
gamma <- c(0.5, -0.3, 0.3)
var.e <- 0.16
beta <- rbind(
  c(0.2, -0.1, 0.1, -0.2),
  c(2.0, -0.1, 0.1, -0.2),
  c(1.0, -0.1, 0.1, -0.2)
)
D <- diag(c(0.16, 0.09, 0.25, 0.16, 0.10, 0.02))
ltriD <- c(0.06, 0.02, 0.04, 0.00, -0.00, 
           0.03, 0.00, -0.06, 0.00,
           0.08, 0.05, 0.01, 
           0.04, -0.00, 
           0.02)
D[lower.tri(D, F)] <- ltriD
D[upper.tri(D, F)] <- t(D)[upper.tri(D, F)]
isSymmetric(D); any(eigen(D)$val<0);det(D)<=0


data <- simData(ntms = 15,
                family = list('gaussian', 'poisson', 'binomial'),
        beta = beta, D = D, var.e = var.e, gamma = gamma, zeta = c(0, -0.2),
        theta = c(-3, 0.2))

source('EM.R')
long.formulas <- list(
  Y.1 ~ time + cont + bin + (1 + time|id),
  Y.2 ~ time + cont + bin + (1 + time|id),
  Y.3 ~ time + cont + bin + (1+time|id)
)
surv.formula <- Surv(survtime, status) ~ bin
EM(long.formulas, surv.formula, data$data, family = list('gaussian', 'poisson', 'binomial'),
   control = list(verbose = T)) -> fit

data <- replicate(100, 
                  simData(ntms = 15,
                          family = list('gaussian', 'poisson', 'binomial'),
                          beta = beta, D = D, var.e = var.e, gamma = gamma, zeta = c(0, -0.2),
                          theta = c(-3, 0.2)),
                  simplify = F)
data <- lapply(data, el)


pb <- utils::txtProgressBar(max=100,style=3)
fits <- vector('list', 100)
for(i in 1:100){
  d <- data[[i]]
  fit <- tryCatch(
    suppressMessages(EM(long.formulas, surv.formula, d, family = list('gaussian', 'poisson', 'binomial'),
                        control = list(verbose = F))),
    error = function(e) NULL
  )
  fits[[i]] <- fit
  utils::setTxtProgressBar(pb, i)
}
do.call(rbind, lapply(fits, function(x) if(!is.null(x)) x$coef$gamma))


# Rustand D ---------------------------------------------------------------
source('simData.R')
D <- matrix(c(0.16, 0.03, 0.02, 0.04, 0.00,
              0.03, 0.09, 0.03, 0.00, -0.06, 
              0.02, 0.03, 0.25, 0.08, 0.05,
              0.04, 0.00, 0.08, 0.16, 0.04, 
              0.00, -0.06, 0.05, 0.04, 0.25),ncol=5,nrow=5)
gamma <- c(0.5, -0.3, 0.3)
var.e <- 0.16
beta <- rbind(
  c(0.2, -0.1, 0.1, -0.2),
  c(2.0, -0.1, 0.1, -0.2),
  c(1.0, -0.1, 0.1, -0.2)
)
data <- simData(ntms = 15,
                family = list('gaussian', 'poisson', 'binomial'),
                beta = beta, D = D, var.e = var.e, gamma = gamma, zeta = c(0, -0.2),
                theta = c(-3, 0.2),
                random.formula = list(~time, ~time, ~1))

source('EM.R')
long.formulas <- list(
  Y.1 ~ time + cont + bin + (1 + time|id),
  Y.2 ~ time + cont + bin + (1 + time|id),
  Y.3 ~ time + cont + bin + (1|id)
)
surv.formula <- Surv(survtime, status) ~ bin
EM(long.formulas, surv.formula, data$data, family = list('gaussian', 'poisson', 'binomial'),
   control = list(verbose = T)) -> fit

# Many fits
data <- replicate(100, 
                  simData(ntms = 15, n = 500,
                          family = list('gaussian', 'poisson', 'binomial'),
                          beta = beta, D = D, var.e = var.e, gamma = gamma, zeta = c(0, -0.2),
                          theta = c(-3, 0.2),
                          random.formula = list(~time, ~time, ~1)),
                  simplify = F)
data <- lapply(data, el)

pb <- utils::txtProgressBar(max=100,style=3)
fits <- vector('list', 100)
for(i in 1:100){
  d <- data[[i]]
  fit <- tryCatch(
    suppressMessages(EM(long.formulas, surv.formula, d, family = list('gaussian', 'poisson', 'binomial'),
                        control = list(verbose = F))),
    error = function(e) NULL
  )
  fits[[i]] <- fit
  utils::setTxtProgressBar(pb, i)
}


# Parse -------------------------------------------------------------------
library(tidyverse)
theme_csda <- function(base_family = "Arial", base_size = 12){
  theme_light(base_family = base_family, base_size = base_size) %+replace%
    theme(
      strip.background = element_blank(),
      strip.text = element_text(size = 11, colour = "black", vjust = 1),
      axis.text = element_text(size = 9),
      axis.title = element_text(size = 11),
      legend.background = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
}
parse.aEM <- function(fit){
  co <- fit$coeffs
  beta <- co$beta
  sigma2 <- co$sigma[[1]]
  gam <- co$gamma
  ze <- co$zeta
  SE <- fit$SE[!grepl('^D', names(fit$SE))]
  Omega <- setNames(c(beta, gam, ze, sigma2), names(SE))
  
  lb <- Omega - qnorm(.975) * SE; ub <- Omega + qnorm(.975) * SE
  return(cbind(Omega, SE, lb, ub))
}
elapsed.aEM <- function(fit, what = 'EM'){
  if(what == 'EM') return(fit$EMtime + fit$postprocess.time)
  else if(what == 'total') return(fit$totaltime)
}

# Tabulation --------------------------------------------------------------
toXdp <- function(x, X) format(round(x, X), nsmall = X)
K <- 3
b <- paste0('beta[', 1:K)
g <- paste0('gamma[', 1:K, ']')
b <- paste0(rep(b, each = 4), 0:3, ']')
nm <- c(b, g, 'zeta', 'sigma2')
targets <- setNames(c(0.2, -0.1, 0.1, -0.2, 2.0, -0.1, 0.1, -0.2, 1.0, -0.1, 0.1, -0.2,
             0.5, -0.3, 0.3, -0.2, 0.16), nm)

all.fits <- do.call(rbind, lapply(fits, function(f){
  ests <- as.data.frame(parse.aEM(f))
  ests$bias <- ests$Omega - targets
  ests$param <- names(targets)
  ests$is.cp <- ests$lb <=  targets & ests$ub >= targets
  ests$target <- targets
  ests
}))

all.fits %>% 
  group_by(param) %>% 
  summarise(.groups = 'keep',
            target = unique(target),
            mean = mean(Omega),
            SE = mean(SE),
            mean.bias = mean(bias),
            SD.bias = sd(bias),
            CP = sum(is.cp)/n(),
            conv = n()/100)

# 3-variate, more individuals? --------------------------------------------
rm(list=ls())
source("EM.R")
D <- matrix(c(0.16, 0.03, 0.02, 0.04, 0.00,
              0.03, 0.09, 0.03, 0.00, -0.06, 
              0.02, 0.03, 0.25, 0.08, 0.05,
              0.04, 0.00, 0.08, 0.16, 0.04, 
              0.00, -0.06, 0.05, 0.04, 0.25),ncol=5,nrow=5)
gamma <- c(0.5, -0.2, 0.3)
beta <- rbind(
  c(0.2, -0.1, 0.1, -0.2),
  c(2.0, -0.1, 0.1, -0.2),
  c(1.0, -1.0, 1.0, -1.0)
)
var.e <- 0.16
data <- simData(ntms = 15, n = 500,
                family = list('gaussian', 'poisson', 'binomial'),
                beta = beta, D = D, var.e = var.e, gamma = gamma, zeta = c(0, -0.2),
                theta = c(-3, 0.2),
                random.formula = list(~time, ~time, ~1))

source('EM.R')
long.formulas <- list(
  Y.1 ~ time + cont + bin + (1 + time|id),
  Y.2 ~ time + cont + bin + (1 + time|id),
  Y.3 ~ time + cont + bin + (1|id)
)
surv.formula <- Surv(survtime, status) ~ bin
EM(long.formulas, surv.formula, data$data, family = list('gaussian', 'poisson', 'binomial'),
   control = list(verbose = T)) -> fit

# 5-variate? --------------------------------------------------------------
rm(list=ls())
gamma <- c(0.5, -0.5, 0.3, -0.3, 0.2)
var.e <- 0.16
beta <- rbind(
  c(0.2, -0.1, 0.1, -0.2),
  -c(0.2, -0.1, 0.1, -0.2),
  c(2.0, -0.1, 0.1, -0.2),
  -c(2.0, -0.1, 0.1, -0.2),
  c(1.0, -0.1, 0.1, -0.2)
)

D <- diag(c(0.25, 0.09, 0.16, 0.04, # Gaussians
            0.25, 0.09, 0.25, 0.05, # Counts
            0.20))                  # Binary
D[2,1] <- D[4,3] <- D[1,2] <- D[4,3] <- 0.03
D[6,5] <- D[8,7] <- D[5,6] <- D[7,8] <- -0.01

isSymmetric(D); any(eigen(D)$val<0);det(D)<=0

source("EM.R")
data <- simData(ntms = 15, n = 500,
                family = list('gaussian', 'gaussian', 'poisson', 'poisson', 'binomial'),
                beta = beta, D = D, var.e = var.e, gamma = gamma, zeta = c(0, -0.2),
                theta = c(-3, 0.2),
                random.formula = list(~time,~time,~time,~time,~1))
data <- data$data
long.formulas <- list(
  Y.1 ~ time + cont + bin + (1 + time|id),
  Y.2 ~ time + cont + bin + (1 + time|id),
  Y.3 ~ time + cont + bin + (1 + time|id),
  Y.4 ~ time + cont + bin + (1 + time|id),
  Y.5 ~ time + cont + bin + (1|id)
)
surv.formula <- Surv(survtime, status) ~ bin
EM(long.formulas, surv.formula, data,
   family = list('gaussian', 'gaussian', 'poisson', 'poisson', 'binomial'),
   control = list(verbose = T, correlated = F, gamma.SE= 'exact')) -> fit
my.summary(fit)

# Many K = 5

# Many fits
data <- replicate(100, 
                  simData(ntms = 15, n = 500,
                          family = list('gaussian', 'gaussian', 'poisson', 'poisson', 'binomial'),
                          beta = beta, D = D, var.e = var.e, gamma = gamma, zeta = c(0, -0.2),
                          theta = c(-3, 0.2),
                          random.formula = list(~time, ~time, ~time, ~time, ~1)),
                  simplify = F)
data <- lapply(data, el)

long.formulas <- list(
  Y.1 ~ time + cont + bin + (1 + time|id),
  Y.2 ~ time + cont + bin + (1 + time|id),
  Y.3 ~ time + cont + bin + (1 + time|id),
  Y.4 ~ time + cont + bin + (1 + time|id),
  Y.5 ~ time + cont + bin + (1|id)
)
surv.formula <- Surv(survtime, status) ~ bin

pb <- utils::txtProgressBar(max=100,style=3)
fits <- vector('list', 100)
for(i in 1:100){
  d <- data[[i]]
  fit <- tryCatch(
    suppressMessages(
      EM(long.formulas, surv.formula, d,
         family = list('gaussian', 'gaussian', 'poisson', 'poisson', 'binomial'),
         control = list(verbose = F, correlated = F, gamma.SE = 'exact'))
    ),
    error = function(e) NULL
  )
  fits[[i]] <- fit
  utils::setTxtProgressBar(pb, i)
}
save(fits, file = '~/Downloads/K5n500-2.RData')



# Parse K=5 ---------------------------------------------------------------
parse.aEM <- function(fit){
  co <- fit$coeffs
  beta <- co$beta
  sigma2 <- c(co$sigma[[1]], co$sigma[[2]])
  gam <- co$gamma
  ze <- co$zeta
  SE <- fit$SE[!grepl('^D', names(fit$SE))]
  Omega <- setNames(c(beta, gam, ze, sigma2), names(SE))
  
  lb <- Omega - qnorm(.975) * SE; ub <- Omega + qnorm(.975) * SE
  return(cbind(Omega, SE, lb, ub))
}

toXdp <- function(x, X) format(round(x, X), nsmall = X)
K <- 5
b <- paste0('beta[', 1:K)
g <- paste0('gamma[', 1:K, ']')
b <- paste0(rep(b, each = 4), 0:3, ']')
nm <- c(b, g, 'zeta', 'sigma2_1', 'sigma2_2')

targets <- setNames(
  c(0.2, -0.1, 0.1, -0.2,
    -0.2, 0.1, -0.1, 0.2,
    2.0, -0.1, 0.1, -0.2,
    -2.0, 0.1, -0.1, 0.2,
    1.0, -0.1, 0.1, -0.2,
    0.5, -0.5, 0.3, -0.3, 0.2,
    -0.2, 0.16, 0.16), nm)

all.fits <- do.call(rbind, lapply(fits, function(f){
  ests <- as.data.frame(parse.aEM(f))
  ests$bias <- ests$Omega - targets
  ests$param <- names(targets)
  ests$is.cp <- ests$lb <=  targets & ests$ub >= targets
  ests$target <- targets
  ests
  if(f$iter != 200) return(ests) else return(0)
}))

N <- sum(unlist(lapply(fits, is.null)))

to.xtab <- all.fits %>% 
  group_by(param) %>% 
  summarise(.groups = 'keep',
            target = unique(target),
            mean = mean(Omega),
            SE = mean(SE),
            mean.bias = mean(bias),
            SD.bias = sd(bias),
            CP = sum(is.cp)/n(),
            conv = n()/99) %>% ungroup 

tab <- to.xtab %>% 
  mutate(
    param = paste0('$\\', param, '=', target,'$')
  ) %>% select(-target) %>% 
  mutate_at(vars(mean, SE, mean.bias, SD.bias), ~ toXdp(.x, 3)) %>% 
  mutate_at('CP', ~ toXdp(.x, 2)) %>% 
  mutate(mean.SD = paste0(mean, ' (', SE, ')'),
         bias.SD = paste0(mean.bias, ' (', SD.bias, ')')) %>% 
  select(param, mean.SD, bias.SD, CP) %>% 
  rename_at(vars(contains('SD')), ~ gsub('\\.SD', ' (SD)', .x))

tab2 <- tab %>% 
  mutate_at('param', ~ gsub('\\[', '_{', .x)) %>% 
  mutate_at('param', ~ gsub('\\]', '}', .x))

xt <- xtable::xtable(tab2)

print(xt,
      include.rownames = FALSE,
      sanitize.text.function = identity)

# ET
quantile(do.call(c, lapply(fits, function(x) if(x$iter<200) x$EMtime + x$postprocess.time)))
