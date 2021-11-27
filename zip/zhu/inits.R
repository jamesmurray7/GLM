#' ####
#' inits.R // Initial conditions for longitudinal and survival sub-models
#'                   for zero-inflated poisson (ZIP) sub-model.
#' --> Same structure as work in Zhu et al. (2018) (int, time, bin)
#' ---
#' Survival part STOPs here.
#' Longitudinal then undergoes MVLME fit
#' ####

if(!'glmmTMB'%in%loadedNamespaces()) library(glmmTMB)
if(!'dplyr'%in%loadedNamespaces()) library(dplyr)

# beta, D inits using lme4 -----------------------------------------

Longit.inits <- function(data){
  fit <- glmmTMB(as.formula(paste0('Y ~ time + bin + (1|id)')),
                 family = poisson, data = data, 
                 ziformula = ~ time + bin + (1|id),
                 control = glmmTMBControl(optCtrl = list(rel.tol = 1e-3)))

  # \beta and D
  beta <- glmmTMB::fixef(fit)$cond
  names(beta) <- paste0('beta', '_', names(beta))

  alpha <- glmmTMB::fixef(fit)$zi
  names(alpha) <- paste0('alpha', '_', names(alpha))
  
  D <- diag(sapply(1:2, function(x) VarCorr(fit)[[x]]$id))
  
  # Checking pos-def. on D, if not then use Matrix::nearPD()$mat
  if(any(eigen(D)$values < 0) || (det(D) <= 0)){
    message("Generated covariance matrix not positive semi-definite")
    message("\n------- -> Transforming... <- -------\n")
    D <- as.matrix(Matrix::nearPD(D, maxit = 1e4)$mat)
  }
  
  list(beta.init = beta,
       alpha.init = alpha,
       D.init = D,
       fit = fit)
}

# Populate REs matrix -----------------------------------------------------
Ranefs <- function(longfits){
  fit <- longfits$fit
  out <- as.matrix(
    cbind(glmmTMB::ranef(fit)$cond$id, glmmTMB::ranef(fit)$zi$id)
  )
  colnames(out) <- c('po_(intercept)', 'zi_(intercept)')
  
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

#ss <- ToStartStop(d)

# Using this StartStop and getting timevarying coxph...
TimeVarCox <- function(data, REs){
  # Prepare data
  ss <- ToStartStop(data)
  REs <- as.data.frame(REs); REs$id <- 1:nrow(REs)
  ss2 <- dplyr::left_join(ss, REs, 'id')
  ss2 <- dplyr::distinct(dplyr::left_join(ss2, data[, c('id', 'bin', 'survtime', 'status')], 'id'))
  
  # Create \gamma_k variables
  gamma <- matrix(0, nrow(ss2), 2)

  gamma[,1] <- ss2[,'po_(intercept)'] 
  gamma[,2] <- ss2[,'zi_(intercept)']
  colnames(gamma) <- c('gamma_1', 'gamma_2')
  
  # And join on ...
  ss3 <- cbind(ss2, gamma)
  # Update this to deal with ties too?
  ss3$status2 <- ifelse(ss3$survtime == ss3$time2, ss3$status, 0)
  
  # Time Varying coxph
  # Reducing only to intercept for this scenario!
  ph <- coxph(Surv(survtime, status) ~ bin + gamma_1 + gamma_2, data = ss3[ss3$time1 == 0, ])
  
  list(inits = coef(ph), l0.init = coxph.detail(ph)$haz, ph = ph)
}


