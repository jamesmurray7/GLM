#' ####
#' inits.R // Initial conditions for longitudinal and survival sub-models
#' ---
#' Survival part STOPs here.
#' Longitudinal then undergoes MVLME fit
#' ####

if(!'glmmTMB'%in%loadedNamespaces()) library(glmmTMB)
if(!'dplyr'%in%loadedNamespaces()) library(dplyr)

# beta, D inits using lme4 -----------------------------------------

Longit.inits <- function(K, data){
  lfitK <- list()
  for(k in 1:K){
    lfitK[[k]] <- glmmTMB(as.formula(paste0('Y.', k, '~ time + cont + bin + (1+time|id)')),
                          family = poisson, data = data, 
                          control = glmmTMBControl(optCtrl = list(rel.tol = 1e-3)))
  }
  
  # \beta and D
  beta <- do.call(c, lapply(1:K, function(k){
    x <- glmmTMB::fixef(lfitK[[k]])$cond
    names(x) <- paste0('beta', k, '_',  names(x))
    x
  }))
  
  D <- as.matrix(Matrix::bdiag(
    lapply(lfitK, function(X){
      matrix(glmmTMB::VarCorr(X)$cond$id, dim(glmmTMB::VarCorr(X)$cond$id))
    })
  ))
  
  # Checking pos-def. on D, if not then use Matrix::nearPD()$mat
  if(any(eigen(D)$values < 0) || (det(D) <= 0)){
    message("Generated covariance matrix not positive semi-definite")
    message("\n------- -> Transforming... <- -------\n")
    D <- as.matrix(Matrix::nearPD(D, maxit = 1e4)$mat)
  }
  
  list(beta.init = beta,
       D.init = D,
       long.fits = lfitK)
}

# Populate REs matrix -----------------------------------------------------
Ranefs <- function(longfits){
  fits <- longfits$long.fits
  K <- length(fits)
  
  # The random effects
  ranefK <- list()
  for(k in 1:K){
    ranefK[[k]] <- as.matrix(glmmTMB::ranef(fits[[k]])$cond$id)
    colnames(ranefK[[k]]) <- paste0(c("intercept_", "slope_"), k)
  }

  REs <- as.data.frame(do.call(cbind, ranefK))
  REs$id <- 1:nrow(REs)
  REs
}

# Survival Inits ----------------------------------------------------------

# Getting data into time1/time2 format...
ToStartStop <- function(data){
  this.subj <- list()
  uids <- unique(data$id)
  
  for(i in uids){
    i.dat <- data[data$id == i, c('id', 'time', 'survtime')]
    df <- as.data.frame(cbind(
      id = i.dat$id,
      time1 = i.dat$time,
      time2 = c(i.dat$time[-1], unique(i.dat$survtime))
    ))
    this.subj[[i]] <- cbind(
      id = i.dat$id,
      time1 = i.dat$time,
      time2 = c(i.dat$time[-1], unique(i.dat$survtime))
    )
  }
  as.data.frame(do.call(rbind, this.subj))
}

#ss <- ToStartStop(d)

# Using this StartStop and getting timevarying coxph...
TimeVarCox <- function(data, REs, fixef.surv = c('cont', 'bin'),
                       survtime = 'survtime', status = 'status'){
  # Prepare data
  ss <- ToStartStop(data)
  ss2 <- dplyr::left_join(ss, REs, 'id')
  ss2 <- dplyr::distinct(dplyr::left_join(ss2, data[, c('id', fixef.surv, survtime, status)], 'id'))
  
  # Create \gamma_k variables
  K <- ncol(REs) %/% 2
  gammaK <- matrix(NA, nrow(ss2), K)
  colnames(gammaK) <- paste0('gamma_', 1:K)
  
  for(k in 1:K){
    ssK <- ss2[, c('time1', paste0(c('intercept_', 'slope_'), k))]
    gammaK[, k] <- ssK[, 2] + ssK[, 1] * ssK[, 3]
  }
  
  # And join on ...
  ss3 <- cbind(ss2, gammaK)
  # Update this to deal with ties too?
  ss3$status2 <- ifelse(ss3$survtime == ss3$time2, ss3$status, 0)
  
  # Time Varying coxph
  # Formula
  timevar.formula <- as.formula(
    paste0('Surv(time1, time2, status2) ~ ', paste0(fixef.surv, collapse = ' + '), ' + ', paste0('gamma_', 1:K, collapse = ' + '))
  )
  ph <- coxph(timevar.formula, data = ss3)
  
  list(inits = coef(ph), l0.init = coxph.detail(ph)$haz, ph = ph)
}

# test <- TimeVarCox(d, REs)





#' #####
#' ARCHIVE
#' #####
# Old version Longit.inits using lme4
# Longit.inits <- function(K, data){
#   lfitK <- list()
#   for(k in 1:K){
#     lfitK[[k]] <- glmer(as.formula(paste0('Y.', k, '~ time + cont + bin + (1+time|id)')),
#                         family = poisson, data = data, 
#                         control = glmerControl(optimizer = 'Nelder_Mead', boundary.tol = 1e-3,
#                                                optCtrl = list(FtolRel = 1e-3, XtolRel = 1e-3)))
#   }
#   
#   # beta and D
#   beta <- do.call(c, lapply(lfitK, lme4::fixef))
#   D <- as.matrix(Matrix::bdiag(
#     lapply(lfitK, function(X){
#       matrix(lme4::VarCorr(X)$id, dim(lme4::VarCorr(X)$id))
#     })
#   ))
#   
#   # Checking pos-def. on D
#   if(any(eigen(D)$values < 0) || (det(D) <= 0)){
#     message("Generated covariance matrix not positive semi-definite")
#     message("\n------- -> Transforming... <- -------\n")
#     D <- as.matrix(Matrix::nearPD(D, maxit = 1e4)$mat)
#   }
#   
#   list(beta.init = beta,
#        D.init = D,
#        long.fits = lfitK)
# }
# Old version Longit.inits using nlme/MASS
# Longit.inits <- function(K, data){
#   lfitK <- list()
#   for(k in 1:K){
#     lfitK[[k]] <- glmmPQL(fixed = as.formula(paste0('Y.', k, '~ time + cont + bin')),
#                           random = ~ 1 + time | id, family = poisson, data = data, niter = 25,
#                           control = lmeControl(opt = 'optim', msTol = 1e-3), verbose = F)
#   }
#   
#   # beta and D
#   beta <- do.call(c, lapply(lfitK, nlme::fixef))
#   D <- as.matrix(Matrix::bdiag(
#     lapply(lfitK, function(X){
#       matrix(nlme::getVarCov(X), dim(nlme::getVarCov(X)))
#     })
#   ))
#   
#   # Checking pos-def. on D
#   if(any(eigen(D)$values < 0) || (det(D) <= 0)){
#     message("Generated covariance matrix not positive semi-definite")
#     message("\n------- -> Transforming... <- -------\n")
#     D <- as.matrix(Matrix::nearPD(D, maxit = 1e4)$mat)
#   }
#   
#   list(beta.init = beta,
#        D.init = D,
#        long.fits = lfitK)
# }