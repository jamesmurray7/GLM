#' ####
#' inits.R // Initial conditions for longitudinal and survival sub-models
#'                   for poisson (univariate currently) sub-model.
#' ####


# Obtain fixed and random effects -----------------------------------------
Longit.inits <- function(long.formula, data, family, dispformula = NULL){
  if(!"formula"%in%class(long.formula)) stop('"long.formula" must be of class "formula"')
  if("function"%in%class(family)) family <- family()$family # step to ensure non-quoted arguments don't throw error.
  
  family <- match.arg(family, c('gaussian', 'binomial', 'poisson', 'negative.binomial'), several.ok = F)
  
  #' Set appropriate family ----
  switch(family, 
         gaussian = f <- gaussian,
         binomial = f <- binomial,
         poisson = f <- poisson,
         negative.binomial = f <- glmmTMB::nbinom2()
  )
  
  # Sanity checks if negative binomial chosen
  if(family == 'negative.binomial' & is.null(dispformula)){ 
    dispformula <- ~ 1
  }else if(family == 'negative.binomial' & !is.null(dispformula)){
    dispformula <- dispformula; warning('Parameterisation of dispersion beyond global intercept not supported.\n')
  }else{
    dispformula <- ~ 1
  }
  
  #' Fit using glmmTMB ----
  fit <- glmmTMB(long.formula,
                 family = f, data = data, dispformula = dispformula,
                 control = glmmTMBControl(optCtrl = list(rel.tol = 1e-3)))

  #' Extract & Return ----
  beta <- glmmTMB::fixef(fit)$cond
  names(beta) <- paste0('beta', '_', names(beta))
  
  D <- glmmTMB::VarCorr(fit)$c$id; dimD <- dim(D)
  D <- matrix(D, nr = dimD[1], nc = dimD[2])
  # Checking pos-def. on D, if not then use Matrix::nearPD()$mat
  if(any(eigen(D)$values < 0) || (det(D) <= 0)){
    message("Generated covariance matrix not positive semi-definite")
    message("\n------- -> Transforming... <- -------\n")
    D <- as.matrix(Matrix::nearPD(D, maxit = 1e4)$mat)
  }
  
  if(family == 'negative.binomial'){
    sigma <- setNames(glmmTMB::sigma(fit), 'theta')
  }else if(family == 'gaussian'){
    sigma <- setNames(glmmTMB::sigma(fit)^2, 'var.e') # This potentially confusing later on!
  }else{
    sigma <- 0
  }
  
  #' Random effects
  b <- as.matrix(glmmTMB::ranef(fit)$cond$id)
    
  list(
    beta.init = beta,
    D.init = D,
    sigma.init = sigma,
    b = b
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
  lhs <- .ToRanefForm(ss3[,'time1'], formulas2$random, q) #sapply(1:q, function(i) ss3[, 'time1']^(i-1))
  gamma <- unname(rowSums(lhs * b[ss3$id,]))
  # gamma2 <- lhs * ss3[,(4:(3+q))]
  # gamma <- unname(rowSums(lhs * rhs))

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

