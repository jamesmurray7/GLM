source('./zip/fresh/_Functions.R')
library(glmmTMB)
# Function to fit, regardless of underlying simulation
fit_mod <- function(data){
  glmmTMB(Y ~ time + cont + bin + (1|id), data = data, family = poisson,
          ziformula = ~ time + (1|id),
          control = glmmTMBControl(optCtrl = list(rel.tol = 1e-3)))
}
# Function to extract parameters
extract_fit <- function(x){
  beta <- fixef(x)$cond
  alpha <- fixef(x)$zi
  D <- sapply(1:2, function(a) VarCorr(x)[[a]]$id)
  out <- c(beta, alpha, D)
  names(out) <- c(paste0('beta', 0:3), paste0('alpha', 0:1), 'D11', 'D22')
  out
}

N <- 10
# Simulations, no survival part -------------------------------------------
beta <- c(1.5, 0.05, 0.33, 0.50)
alpha <- c(-0.5, 0.25)
D <- diag(c(.5^2, .15^2))
n <- 250
ntms <- 15
data.nosurv <- replicate(N, 
                  simData_zip(n, ntms, beta, alpha, D), simplify = F)

# Simulations, survival part ----------------------------------------------
beta <- c(1.5, 0.05, 0.33, 0.50)
alpha <- c(-0.5, 0.25)
D <- diag(c(.5^2, .15^2))
gamma <- .5
n <- 250
ntms <- 15
data.surv.L <- replicate(N, 
                         simData_zip_joint(n, ntms, beta, alpha, D,
                                           theta = c(-6, 0.2), gamma = gamma), simplify = F)
data.surv.H <- replicate(N, 
                           simData_zip_joint(n, ntms, beta, alpha, D,
                                             theta = c(-4.5, 0.2), gamma = gamma), simplify = F)

# Fit and compare ---------------------------------------------------------

fits.nosurv <- fits.surv.H <- fits.surv.L <- list()
pb <- utils::txtProgressBar(max = N, style = 3)
for(m in 1:N){
  fits.nosurv[[m]] <- tryCatch(fit_mod(data.nosurv[[m]]$data), error = function(e) NULL)
  fits.surv.H[[m]] <- tryCatch(fit_mod(data.surv.H[[m]]$data), error = function(e) NULL)
  fits.surv.L[[m]] <- tryCatch(fit_mod(data.surv.L[[m]]$data), error = function(e) NULL)
  utils::setTxtProgressBar(pb, m)
}

ll <- function(a) sum(unlist(lapply(a, is.null))) # All appear to have fit ...
ll(fits.nosurv); ll(fits.surv.H); ll(fits.surv.L)


# Plot --------------------------------------------------------------------
library(tidyverse)
df.nosurv <- data.frame(a = 'No survival part', do.call(rbind, lapply(fits.nosurv, extract_fit)))
df.surv.H <- data.frame(a = '50% failure surv', do.call(rbind, lapply(fits.surv.H, extract_fit)))
df.surv.L <- data.frame(a = '30% failure surv', do.call(rbind, lapply(fits.surv.L, extract_fit)))

targets <- data.frame(name = names(df.nosurv)[-1],
                      target = c(beta, alpha, diag(D)))

tdf <- rbind(df.nosurv, df.surv.H, df.surv.L)
df %>% 
  pivot_longer(-a) %>%
  left_join(., targets , 'name') %>% 
  ggplot(aes(x = value, colour = a)) +
  geom_vline(aes(xintercept = target), lty = 3) + 
  geom_density() + 
  facet_wrap(~name, scales = 'free')




