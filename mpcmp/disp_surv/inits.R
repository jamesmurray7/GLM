#' ####
#' inits.R // Initial conditions for longitudinal and survival sub-models
#'                   for poisson (univariate currently) sub-model.
#' ####


# Obtain fixed and random effects -----------------------------------------
Longit.inits <- function(long.formula, data){
 
  
  #' Fit using glmmTMB ----
  fit <- glmmTMB(long.formula,
                 family = poisson, data = data)

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
  
  #' Random effects
  b <- as.matrix(glmmTMB::ranef(fit)$cond$id)
  
  list(
    beta.init = beta,
    D.init = D,
    sigma.init = sigma,
    b = b                  
  )
}

# Dispersion --------------------------------------------------------------
source('disp_inits.R')
get.delta.inits <- function(dmats, beta, b, method, summax = NULL, verbose = F, min.profile.length = 1, max.val = Inf){
  
  #' Data objects ----
  X <- dmats$X; Y <- dmats$Y; Z <- dmats$Z # Longitudinal data matrices
  G <- dmats$G                             # Dispersion data matrix
  N <- length(b)
  
  if(is.null(summax)) summax <- ceiling(max(sapply(Y, max)) * 2) else summax <- summax
  mus <- mapply(function(X, Z, b) exp(X %*% beta + Z %*% b), X = X, Z = Z, b = b)
  
  if(verbose) message('Obtaining initial estimate for dispersion parameter delta...\n')
  a <- proc.time()[3]
  if(method == 'bobyqa'){
    raw <- find.deltas.bobyqa(Y,mus, summax, verbose, min.profile.length)
  }else{
    raw <- find.deltas.optim(Y, mus, summax, verbose, min.profile.length, max.val)  
  }
  b <- proc.time()[3]
  
  # Adding to discount instances where estimate doesn't at all move.
  o <- raw
  o[round(o, 3) < -max.val] <- -max.val; o[round(o, 3) > max.val] <- max.val
  o[is.na(o)] <- 0.
  
  # Return
  list(
    subject.estimates = raw,
    truncated.estimates = o,
    median.estimate = median(raw, na.rm = T),
    mean.estimate = mean(raw, na.rm = T),
    median.cut.estimate = median(o),
    mean.cut.estimate = mean(o),
    IQR.estimates = IQR(raw, na.rm = T),
    sd.estimates = sd(raw, na.rm = T),
    time = round(b - a,3)
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

.ToDispForm <- function(time, data, disp.formula){
  time <- as.data.frame(time)
  names(time)[2] <- 'time'
  df <- merge(time, 
              data[,names(data)!='time'], 
              'id')
  df <- df[!duplicated(df[,c('id', 'time')]),]
  model.matrix(disp.formula, df)
}
 
# Using this StartStop and getting timevarying coxph...
TimeVarCox <- function(data, b, ph, formulas, disp.formula, deltas, max.val){
  # Prepare data
  ss <- .ToStartStop(data); q <- ncol(b) # send to Start-Stop (ss) format
  REs <- as.data.frame(b); REs$id <- 1:nrow(b)
  ss2 <- merge(ss, REs, 'id')
  ss3 <- merge(ss2, data[, c('id', colnames(ph$x), 'survtime', 'status')], 'id')
  ss3 <- ss3[!duplicated.matrix(ss3), ]
  
  # Create gamma variable for SRE part...
  lhs <- .ToRanefForm(ss3[,'time1'], formulas$random, q) #sapply(1:q, function(i) ss3[, 'time1']^(i-1))
  gamma1 <- unname(rowSums(lhs * b[ss3$id,]))
  # Create gamma variable for dispersion part...
  lhs <- .ToDispForm(ss3[,c('id', 'time1')], data, disp.formula)
  print(head(lhs))
  dd <- do.call(rbind, deltas)
  dd <- apply(dd, 2, function(d) ifelse(abs(d) >= max.val, 0, d)) # If beyond truncation amount then don't contribute.
  gamma2 <- unname(rowSums(lhs * dd[ss3$id,]))
  
  # And join on ...
  ss3 <- cbind(ss3, gamma1, gamma2)
  ss3$status2 <- ifelse(ss3$survtime == ss3$time2, ss3$status, 0)
  
  # Time Varying coxph
  # Formula
  timevar.formula <- as.formula(
    paste0('Surv(time1, time2, status2) ~ ', paste0(colnames(ph$x), collapse = ' + '), ' + gamma1 + gamma2')
  )
  ph <- coxph(timevar.formula, data = ss3)
  
  list(inits = coef(ph), l0.init = coxph.detail(ph)$haz, ph = ph)
}

