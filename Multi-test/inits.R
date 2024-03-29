#' ####
#' inits.R // Initial conditions for longitudinal and survival sub-models.
#' ####


# Obtain fixed and random effects -----------------------------------------
Longit.inits <- function(long.formula, data, family, dispformula = NULL){
  lapply(long.formula, function(x) if(!"formula"%in%class(x)) stop('"long.formula" must be of class "formula"'))
  family.form <- lapply(family, function(f){
    
    if("function"%in%class(f)) f <- f()$family # step to ensure non-quoted arguments don't throw error.
    ff <- match.arg(f, c('gaussian', 'binomial', 'poisson', 'negative.binomial'), several.ok = F)
    
    #' Set appropriate family ----
    switch(ff, 
           gaussian = family <- gaussian,
           binomial = family <- binomial,
           poisson = family <- poisson,
           negative.binomial = family <- glmmTMB::nbinom2()
    )
    family
  })
  
  K <- length(long.formula)
  if(K!=length(family)) stop('Uneven long.formula and family supplied.')
  #' Fit using glmmTMB ----
  fits <- lapply(1:K, function(k){
    dispformula <- ~ 1
    fit <- glmmTMB(long.formula[[k]],
                   family = family.form[[k]], data = data, dispformula = dispformula,
                   control = glmmTMBControl(optCtrl = list(rel.tol = 1e-3)))
    fit
  })
  
  #' Extract -----                (there is likely a way using mapply to do this all-in-one!)
  #' Biomarker names
  markers <- unlist(lapply(1:K, function(k) names(fits[[k]]$modelInfo$respCol)))
  #' Fixed-effects
  beta <- do.call(c, lapply(1:K, function(k){
    coef <- glmmTMB::fixef(fits[[k]])$cond
    names(coef) <- paste0(markers[k], '_', names(coef))
    coef
  }))
  #' Random effects
  Ds <- lapply(1:K, function(k){
    D <- glmmTMB::VarCorr(fits[[k]])$c$id; dimD <- dim(D)
    D <- matrix(D, nr = dimD[1], nc = dimD[2])
    # Checking pos-def. on this D, if not then use Matrix::nearPD()$mat
    if(any(eigen(D)$values < 0) || (det(D) <= 0)){
      message("Generated covariance matrix not positive semi-definite, occurred for ", markers[k], ".")
      message("\n------- -> Transforming... <- -------\n")
      D <- as.matrix(Matrix::nearPD(D, maxit = 1e4)$mat)
    }
    D
  })
  D <- as.matrix(Matrix::bdiag(Ds))
  
  #' Dispersions
  sigma <- lapply(1:K, function(k){
    if("function"%in%class(family[[k]])) f <- family[[k]]()$family else f <- family[[k]]
    f <- match.arg(f, c('gaussian', 'binomial', 'poisson', 'negative.binomial'), several.ok = F)
    if(f=='negative.binomial'|f=='nbinom2'){
      out <- setNames(glmmTMB::sigma(fits[[k]]), 'theta')
    }else if(f=='gaussian'){
      out <- setNames(glmmTMB::sigma(fits[[k]])^2, 'var.e')
    }else{
      out <- 0
    }
    out
  })
  
  #' Random effects
  b <- do.call(cbind, lapply(1:K, function(k){
    b <- as.matrix(glmmTMB::ranef(fits[[k]])$cond$id)
    colnames(b) <- paste0(markers[k], '_', colnames(b))
    b
  }))
    
  
  #' Populate off-block-diagonal terms in D 
  off.inds <- which(D == 0, arr.ind = T)
  D[lower.tri(D, F)] <- cov(b)[lower.tri(cov(b), F)]
  D[upper.tri(D, F)] <- t(D)[upper.tri(D, F)]
    
  list(
    beta.init = beta,
    D.init = D,
    sigma.init = sigma,
    b = b,
    responses = markers,
    off.inds = off.inds
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

.ToRanefForm <- function(time, random.formula){
  if(attr(random.formula, 'special') == 'none'){
    out <- model.matrix(as.formula(paste0('~', random.formula)), as.data.frame(time))
  }else if(attr(random.formula, 'special') == 'spline'){
    out <- model.matrix(as.formula(paste0('~', random.formula)), as.data.frame(time))
  }
  as.data.frame(out)
}

# Using this StartStop and getting timevarying coxph...
TimeVarCox <- function(data, b, ph, formulas, b.inds){
  # Prepare data
  ss <- .ToStartStop(data); q <- ncol(b) # send to Start-Stop (ss) format
  REs <- as.data.frame(b); REs$id <- 1:nrow(b); K <- length(b.inds)
  ss2 <- merge(ss, REs, 'id')
  ss3 <- merge(ss2, data[, c('id', colnames(ph$x), 'survtime', 'status')], 'id')
  ss3 <- ss3[!duplicated.matrix(ss3), ]
  
  # Create gamma variable
  lhs <- lapply(formulas, function(x) .ToRanefForm(ss3[,'time1'], x$random))
  gamma <- lapply(1:K, function(k){
    unname(rowSums(lhs[[k]] * b[ss3$id, b.inds[[k]]]))
  })
  gamma <- do.call(cbind, gamma)
  colnames(gamma) <- paste0('gamma_', 1:K)
  
  # And join on ...
  ss3 <- cbind(ss3, gamma)
  # Update this to deal with ties too?
  ss3$status2 <- ifelse(ss3$survtime == ss3$time2, ss3$status, 0)
  
  # Time Varying coxph
  # Formula
  timevar.formula <- as.formula(
    paste0('Surv(time1, time2, status2) ~ ', paste0(colnames(ph$x), collapse = ' + '), ' + ', paste0('gamma_', 1:K, collapse = ' + '))
  )
  ph <- coxph(timevar.formula, data = ss3)
  
  list(inits = coef(ph), l0.init = coxph.detail(ph)$haz, ph = ph)
}

