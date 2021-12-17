load('~/Downloads/mix-sim-data.RData')
data1 <- data[[1]]

library(JMbayes2)
fitn <- function(x){
  ph <- coxph(Surv(survtime, status) ~ cont + bin, data = x$surv.data)
  m1 <- lme(Y.1 ~ time + cont + bin, random = ~ time |id, data = x$data)
  m2 <- mixed_model(Y.2 ~ time + cont + bin, random = ~ time | id, data = x$data, family = binomial())
  m3 <- mixed_model(Y.3 ~ time + cont + bin, random = ~ time | id, data = x$data, family = poisson())
  
  M <- list(m1, m2, m3)
  fit <- tryCatch(
    jm(ph, M, time_var = 'time', 
        n_iter = 11000, n_burnin = 1000L, n_thin = 5L, n_chains = 1L, cores = 1),
    error = function(e) NULL
    )
  return(fit)
}
vech <- function(x) x[lower.tri(x, T)]
jmb.extract <- function(jmb.fit){
  sj <- summary(jmb.fit)
  G <- sj$Outcome1
  B <- sj$Outcome2
  P <- sj$Outcome3
  surv <- sj$Surv
  
  # Setting names to match my fit for ease of comparison.
  beta_rhs <- c('_(Intercept)', '_time', '_cont', '_bin')
  # 1. Gaussian
  G_beta <- G[1:4,1]; names(G_beta) <- paste0('G', beta_rhs)
  var.e <- G[5,1]^2; names(var.e) <- 'var.e'
  # CI
  G_beta95 <- G[1:4, 3:4]; row.names(G_beta95) <- names(G_beta)
  var.e95 <- G[5, 3:4]^2; row.names(var.e95) <- 'var.e'
  # 2. Binomial
  B_beta <- B[1:4, 1]; names(B_beta) <- paste0('B', beta_rhs)
  # CI
  B_beta95 <- B[1:4, 3:4];  row.names(B_beta95) <- names(B_beta)
  # 3. Poisson
  P_beta <- P[1:4, 1]; names(P_beta) <- paste0('P', beta_rhs)
  # CI
  P_beta95 <- P[1:4, 3:4];  row.names(P_beta95) <- names(P_beta)
  # 4. Survival
  eta <- surv[1:2, 1]; names(eta) <- c('surv_cont', 'surv_bin')
  gamma <- surv[3:5, 1]; names(gamma) <- paste0('gamma_', 1:3)
  # CI
  eta95 <- surv[1:2, 3:4]; row.names(eta95) <- names(eta)
  gamma95 <- surv[3:5, 3:4]; row.names(gamma95) <- names(gamma)
  # 5. Covariance D
  D <- vech(sj$D)
  names(D) <- paste0('[', apply(which(lower.tri(sj$D, T), arr.ind = T), 1, paste0, collapse = ', '), ']')
  
  # Pack and return
  estimates <- list(
    D,
    c(G_beta, B_beta, P_beta),var.e,
    gamma, eta
  )
  CIs <- list(
    rbind(G_beta95, B_beta95, P_beta95), var.e95,
    gamma95, eta95
  )
  return(list(estimates, CIs, jmb.fit$running_time[3]))
}

pb <- utils::txtProgressBar(max = 100, style = 3)
fit <- list()
for(i in 1:100){
  fit[[i]] <- fitn(data1[[i]])
  utils::setTxtProgressBar(pb, i)
}

# Deparse results

jmb.fits <- lapply(fit, jmb.extract)
save(jmb.fits, file = '~/Downloads/jmbfits.RData')
