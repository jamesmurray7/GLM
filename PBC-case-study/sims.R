rm(list=ls(
  
))
source('EM.R')

long.formula <- Y ~ time + cont + bin + (1 + time|id)
surv.formula <- Surv(survtime, status) ~ cont + bin


M <- 100
G <- replicate(M, simData(family = gaussian, theta = c(-2.5, .1)), simplify = F)
B <- replicate(M, simData(family = binomial, theta = c(-2.5, .1)), simplify = F)
P <- replicate(M, simData(family = poisson, theta = c(-2.5, .1)), simplify = F)
N <- replicate(M, simData(family = "negative.binomial", theta = c(-2.5, .1), disp = 2), simplify = F)

pb <- utils::txtProgressBar(max = M, style = 3)
Gfit <- Bfit <- Pfit <- Nfit <- vector('list', M)
for(i in 1:M){
  Gfit[[i]] <- tryCatch(suppressMessages(EM(long.formula, surv.formula, G[[i]]$data, gaussian)), error = function(e) NULL)
  Bfit[[i]] <- tryCatch(suppressMessages(EM(long.formula, surv.formula, B[[i]]$data, binomial)), error = function(e) NULL)
  Pfit[[i]] <- tryCatch(suppressMessages(EM(long.formula, surv.formula, P[[i]]$data, poisson)), error = function(e) NULL)
  Nfit[[i]] <- tryCatch(suppressMessages(EM(long.formula, surv.formula, N[[i]]$data, "negative.binomial")), error = function(e) NULL)
  
  utils::setTxtProgressBar(pb, i)
}

fits <- list(`G` = Gfit, `B` = Bfit, `P` = Pfit, `N` = Nfit)
save(fits, file = '~/Downloads/allfits.RData')
