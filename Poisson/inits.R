#' ####
#' inits.R // Initial conditions for longitudinal and survival sub-models
#'                   for poisson (univariate currently) sub-model.
#' ---
#' Survival part STOPs here.
#' Longitudinal then undergoes MVLME fit
#' ####

if(!'glmmTMB'%in%loadedNamespaces()) library(glmmTMB)
if(!'dplyr'%in%loadedNamespaces()) library(dplyr)

# beta, D inits using lme4 -----------------------------------------

Longit.inits <- function(data){
  fit <- glmmTMB(as.formula(paste0('Y ~ time + cont + bin + (1 + time|id)')),
                 family = poisson, data = data, dispformula = ~ 1,
                 control = glmmTMBControl(optCtrl = list(rel.tol = 1e-3)))

  # \beta and D
  beta <- glmmTMB::fixef(fit)$cond
  names(beta) <- paste0('beta', '_', names(beta))
  
  D <- matrix(glmmTMB::VarCorr(fit)$c$id, 2, 2)
  
  # Checking pos-def. on D, if not then use Matrix::nearPD()$mat
  if(any(eigen(D)$values < 0) || (det(D) <= 0)){
    message("Generated covariance matrix not positive semi-definite")
    message("\n------- -> Transforming... <- -------\n")
    D <- as.matrix(Matrix::nearPD(D, maxit = 1e4)$mat)
  }
  
  list(beta.init = beta,
       D.init = D,
       theta.init = glmmTMB::sigma(fit),
       fit = fit)
}

# Populate REs matrix -----------------------------------------------------
Ranefs <- function(longfits){
  fit <- longfits$fit
  out <- as.matrix(glmmTMB::ranef(fit)$cond$id)
  colnames(out) <- c('(intercept)', '(slope)')
  out
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
  REs <- as.data.frame(REs); REs$id <- 1:nrow(REs)
  ss2 <- dplyr::left_join(ss, REs, 'id')
  ss2 <- dplyr::distinct(dplyr::left_join(ss2, data[, c('id', 'cont', 'bin', 'survtime', 'status')], 'id'))
  
  # Create \gamma_k variables
  gamma <- matrix(NA, nrow(ss2), 1)
  gamma[,1] <- ss2[,'(intercept)'] + ss2[,'(slope)'] * ss2[, 'time1']
  
  # And join on ...
  ss3 <- cbind(ss2, gamma)
  # Update this to deal with ties too?
  ss3$status2 <- ifelse(ss3$survtime == ss3$time2, ss3$status, 0)
  
  # Time Varying coxph
  # Formula
  timevar.formula <- as.formula(
    paste0('Surv(time1, time2, status2) ~ ', paste0(c('cont', 'bin'), collapse = ' + '), ' + gamma')
  )
  ph <- coxph(timevar.formula, data = ss3)
  
  list(inits = coef(ph), l0.init = coxph.detail(ph)$haz, ph = ph)
}

