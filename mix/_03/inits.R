#' ####
#' inits.R // Initial conditions for longitudinal and survival sub-models
#'    for a mixture of guassian (1); binomial (2) and poisson(3) sub-models
#' ---
#' Survival part STOPs here.
#' Longitudinal then undergoes MVLME fit
#' ---
#' This directory an implementation of Rustand et al. Simulation (ArXiv 220306256).
#' ####

if(!'glmmTMB'%in%loadedNamespaces()) library(glmmTMB)
if(!'dplyr'%in%loadedNamespaces()) library(dplyr)

# beta, D inits using lme4 -----------------------------------------
Longit.inits <- function(data){
  fits <- pf <- list()
  for(i in 1:3){
    if(i == 1){f = gaussian; prefix = 'G'; RE <- paste0('(1 + time|id)')}
    if(i == 2){f = binomial; prefix = 'B'; RE <- paste0('(1|id)')}
    if(i == 3){f = poisson;  prefix = 'P'; RE <- paste0('(1 + time|id)')}
    
    fits[[i]] <- glmmTMB(as.formula(paste0('c(Y.', i, ') ~ time + cont + bin +', RE)),
                         family = f, data = data,
                         control = glmmTMBControl(optCtrl = list(rel.tol = 1e-3)))
    pf[[i]] <- prefix
  }
  
  # Extract beta
  beta <- do.call(c, lapply(1:3, function(i){
    x <- glmmTMB::fixef(fits[[i]])$cond
    names(x) <- paste0(pf[[i]], '_', names(x))
    x
  }))
  
  # Extract sigma^2
  var.e <- sigma(fits[[1]])^2
  names(var.e) <- 'var.e'
  
  # Extract D
  Ds <- lapply(fits, glmmTMB::VarCorr)
  D <- as.matrix(Matrix::bdiag(lapply(1:3, function(i){
    matrix(Ds[[i]]$cond$id, dim(Ds[[i]]$cond$id))
  })))
  
  # Checking pos-def. on D, if not then use Matrix::nearPD()$mat
  if(any(eigen(D)$values < 0) || (det(D) <= 0)){
    message("Generated covariance matrix not positive semi-definite")
    message("\n------- -> Transforming... <- -------\n")
    D <- as.matrix(Matrix::nearPD(D, maxit = 1e4)$mat)
  }
  
  list(beta.init = beta,
       var.e.init = var.e,
       D.init = D,
       fits = fits)
}

# Populate REs matrix -----------------------------------------------------
extractREformula <- function(fit){
  one <- el(strsplit(paste0(fit$call$formula)[3], '\\('), 1)[2]
  two <- el(strsplit(one, '\\|'), 1)[1]
  three <- unname(sapply(el(strsplit(two, '\\+'), 1), trimws))
  gsub('1', '(Intercept)', three)
} # non-lazy fix to different RE structures.

Ranefs <- function(longfits){
  pf <- unique(substr(names(longfits$beta.init), 1, 2))
  fits <- longfits$fits
  out <- lapply(seq(length(pf)), function(i){
    x <- as.matrix(glmmTMB::ranef(fits[[i]])$cond$id)
    colnames(x) <- paste0(pf[i], extractREformula(fits[[i]]))
    x
  })
  as.matrix(do.call(cbind, out))
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

# Using this StartStop and getting timevarying coxph...
TimeVarCox <- function(data, REs){
  # Prepare data
  ss <- ToStartStop(data)
  pfs <- unique(substr(colnames(REs), 1, 2))
  REs <- as.data.frame(REs); REs$id <- 1:nrow(REs)
  ss2 <- dplyr::left_join(ss, REs, 'id')
  ss2 <- dplyr::distinct(dplyr::left_join(ss2, data[, c('id', 'survtime', 'status')], 'id'))
  
  # Create \gamma_k variables
  
  gamma <- matrix(NA, nrow(ss2), length(pfs))
  for(i in seq_along(pfs)){
    if (i != 2){
      gamma[, i] <- ss2[, paste0(pfs[i], '(Intercept)')] + ss2[,paste0(pfs[i], 'time')] * ss2[, 'time1']
    }else{
      gamma[, i] <- ss2[, paste0(pfs[i], '(Intercept)')] # hardcode binomial part
    }
  }
  colnames(gamma) <- paste0('gamma_', 1:length(pfs))
  
  # And join on ...
  ss3 <- cbind(ss2, gamma)
  # Update this to deal with ties too?
  ss3$status2 <- ifelse(ss3$survtime == ss3$time2, ss3$status, 0)
  
  # Time Varying coxph
  # Formula
  timevar.formula <- as.formula(
    paste0('Surv(time1, time2, status2) ~ ', paste0('gamma_', 1:length(pfs), collapse = ' + '))
  )
  ph <- coxph(timevar.formula, data = ss3)
  
  list(inits = coef(ph), l0.init = coxph.detail(ph)$haz, ph = ph)
}

