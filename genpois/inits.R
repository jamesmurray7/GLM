#' ####
#' inits.R // Initial conditions for longitudinal and survival sub-models
#'                   for poisson sub-model.
#' ####


# Obtain fixed and random effects -----------------------------------------
Longit.inits <- function(long.formula, disp.formula, data, genpois.inits){
 
  
  #' Fit using glmmTMB ----
  if(!genpois.inits){
    fit <- glmmTMB(long.formula,
                   family = poisson, data = data)
    
    #' Extract & Return ----
    beta <- glmmTMB::fixef(fit)$cond
    names(beta) <- paste0('beta', '_', names(beta))
    
    phi.init <- NULL
  }else{ # Else fit with generalised poisson
    fit <- suppressWarnings(glmmTMB(long.formula, data, family = genpois,
                   dispformula = disp.formula))
    #' Extract parameters
    beta <- glmmTMB::fixef(fit)$cond
    names(beta) <- paste0('beta', '_', names(beta))
    
    #' Dispersion is characterised as -{cmp disp}
    phi.init <- exp(glmmTMB::fixef(fit)$disp/2) - 1
    phi.init <- setNames(phi.init, paste0('phi_', names(phi.init)))
  }
  
  D <- glmmTMB::VarCorr(fit)$c$id; dimD <- dim(D)
  D <- matrix(D, nr = dimD[1], nc = dimD[2])
  # Checking pos-def. on D, if not then use Matrix::nearPD()$mat
  if(any(eigen(D)$values < 0) || (det(D) <= 0)){
    message("Generated covariance matrix not positive semi-definite")
    message("\n------- -> Transforming... <- -------\n")
    D <- as.matrix(Matrix::nearPD(D, maxit = 1e4)$mat)
  }
  
  #' Random effects
  b <- as.matrix(glmmTMB::ranef(fit)$cond$id)
  
  list(
    beta.init = beta,
    D.init = D,
    sigma.init = sigma,
    b = b,
    phi.init = phi.init
  )
}

# Survival Inits ----------------------------------------------------------
# Getting data into time1/time2 format...
.ToStartStop <- function(data){
  this.subj <- list()
  uids <- unique(data$id)
  for(i in uids){
    i.dat <- data[data$id == i, c('id', 'time', 'survtime')]
    this.subj[[i]] <- cbind(
      id = i.dat$id,
      time1 = i.dat$time,
      time2 = c(i.dat$time[-1], unique(i.dat$survtime))
    )
  }
  as.data.frame(do.call(rbind, this.subj))
}

.ToRanefForm <- function(time, random.formula, q){
  if(attr(random.formula, 'special') == 'none'){
    out <- sapply(1:q, function(i) time^(i - 1))
  }else if(attr(random.formula, 'special') == 'spline'){
    out <- model.matrix(as.formula(paste0('~', random.formula)), as.data.frame(time))
  }
  as.data.frame(out)
}
 
# Using this StartStop and getting timevarying coxph...
TimeVarCox <- function(data, b, ph, formulas){
  # Prepare data
  ss <- .ToStartStop(data); q <- ncol(b) # send to Start-Stop (ss) format
  REs <- as.data.frame(b); REs$id <- 1:nrow(b)
  ss2 <- merge(ss, REs, 'id')
  ss3 <- merge(ss2, data[, c('id', colnames(ph$x), 'survtime', 'status')], 'id')
  ss3 <- ss3[!duplicated.matrix(ss3), ]
  
  # Create gamma variable
  lhs <- .ToRanefForm(ss3[,'time1'], formulas$random, q) 
  gamma <- unname(rowSums(lhs * b[ss3$id,]))
  # And join on ...
  ss3 <- cbind(ss3, gamma)
  # Update this to deal with ties too?
  ss3$status2 <- ifelse(ss3$survtime == ss3$time2, ss3$status, 0)
  
  # Time Varying coxph
  # Formula
  timevar.formula <- as.formula(
    paste0('Surv(time1, time2, status2) ~ ', paste0(colnames(ph$x), collapse = ' + '), ' + gamma')
  )
  ph <- coxph(timevar.formula, data = ss3)
  
  list(inits = coef(ph), l0.init = coxph.detail(ph)$haz, ph = ph)
}

