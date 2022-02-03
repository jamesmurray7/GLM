#' ####
#' inits.R // Initial conditions for longitudinal and survival sub-models.
#' ---
#' Survival part STOPs here.
#' Longitudinal then undergoes MVLME fit
#' ####

library(nlme)
library(dplyr)

# var.e, beta, D inits using lme4 -----------------------------------------
Longit.inits <- function(K, data, degree, basis = 'survival'){
  lfitK <- list()
  for(k in 1:K){
    if(basis == 'survival'){
      fixed.formula <- as.formula(paste0('Y.', k, ' ~ ', paste0('basis', 1:degree, collapse = ' + '), ' + cont + bin'))
      random.formula <- as.formula(paste0(' ~ ', paste0('basis', 1:degree, collapse = ' + '), ' | id'))
    }else if(basis == 'auto'){
      fixed.formula <- as.formula(paste0('Y.', k, ' ~ bs(time, degree = ', degree, ') + cont + bin'))
      random.formula <- as.formula(paste0(' ~ bs(time, degree = ', degree, ')|id'))
    }else{
      stop('Basis used to obtain univariate model fits must be one of c("survival", "auto").')
    }
    lfitK[[k]] <- lme(fixed = fixed.formula,
                      random = random.formula, data = data,
                      method = "ML",
                      control = nlme::lmeControl(opt = "optim", msTol = 1e-3))
  }
  
  # The three items
  var.e <- do.call(c, lapply(1:K, function(k){
    x <- lfitK[[k]]$sigma
    names(x) <- paste0('var.e_', k)
    x
  }))^2
  
  beta <- do.call(c, lapply(1:K, function(k){
    x <- nlme::fixef(lfitK[[k]])
    names(x) <- paste0('beta', k, '_', names(x))
    x
  }))
  
  # D 
  D <- as.matrix(Matrix::bdiag(
    lapply(lfitK, function(X){
      matrix(getVarCov(X), dim(getVarCov(X)))
    })
  ))
  
  if(any(eigen(D)$values < 0) || (det(D) <= 0)){
    message("Generated covariance matrix not positive semi-definite")
    message("\n------- -> Transforming... <- -------\n")
    D <- Matrix::nearPD(D, maxit = 1e4)$mat
  }
  
  list(var.e.init = var.e, 
       beta.init = beta,
       D.init = D,
       long.fits = lfitK)
}
# inits.long <- Longit.inits(3, 3, data)

# Populate REs matrix -----------------------------------------------------
Ranefs <- function(longfits){
  fits <- longfits$long.fits
  K <- length(fits)
  
  # The random effects
  ranefK <- list()
  for(k in 1:K){
    re <- ranef(fits[[k]])
    ranefK[[k]] <- as.matrix(re)
    colnames(ranefK[[k]]) <- paste0('b', k, '_', colnames(re))
  }
  # collate and return
  as.matrix(do.call(cbind, ranefK))
}

# Survival Inits ----------------------------------------------------------


# Maybe deprecated? Updating to use basis below... ------------------------

# Getting data into time1/time2 format...
ToStartStop <- function(data){
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

#ss <- ToStartStop(d)

# Using this StartStop and getting timevarying coxph...
TimeVarCox <- function(data, REs, nK = 3, fixef.surv = c('cont', 'bin'),
                       survtime = 'survtime', status = 'status'){
  # Prepare data
  ss <- ToStartStop(data)
  REs <- as.data.frame(REs); REs$id <- 1:nrow(REs)
  ss2 <- dplyr::left_join(ss, REs, 'id')
  ss2 <- dplyr::distinct(dplyr::left_join(ss2, data[, c('id', fixef.surv, survtime, status)], 'id'))
  
  # Create \gamma_k variables
  gammaK <- matrix(NA, nrow(ss2), nK)
  colnames(gammaK) <- paste0('gamma_', 1:nK)
  
  for(k in 1:nK){
    ii <- grepl(paste0('time1|^b', k), x = colnames(ss2))
    ssK <- as.matrix(ss2[, ii])
    time.ii <- grepl('time1', colnames(ssK))
    inte.ii <- grepl('\\(Intercept\\)', colnames(ssK))
    gammaK[, k] <- ssK[, inte.ii] + ssK[, time.ii, drop = F] * apply(ssK[, !time.ii & !inte.ii], 1, sum)
  }
  
  # And join on ...
  ss3 <- cbind(ss2, gammaK)
  # Update this to deal with ties too?
  ss3$status2 <- ifelse(ss3$survtime == ss3$time2, ss3$status, 0)
  
  # Time Varying coxph
  # Formula
  timevar.formula <- as.formula(
    paste0('Surv(time1, time2, status2) ~ ', paste0(fixef.surv, collapse = ' + '), ' + ', paste0('gamma_', 1:nK, collapse = ' + '))
  )
  ph <- coxph(timevar.formula, data = ss3)
  
  list(inits = coef(ph), l0.init = coxph.detail(ph)$haz, ph = ph)
}


# Version using spline bases ----------------------------------------------

# Getting data into time1/time2 format...
ToStartStop2 <- function(data, basis){
  this.subj <- list()
  uids <- unique(data$id)
  deg <- attr(basis, 'degree')
  for(i in uids){
    i.dat <- data[data$id == i, c('id', 'time', 'survtime')]
    this.subj[[i]] <- cbind(
      id = i.dat$id,
      time1 = i.dat$time,
      time2 = c(i.dat$time[-1], unique(i.dat$survtime))
    )
  }
  tss <- as.data.frame(do.call(rbind, this.subj))
  suppressWarnings(structure(as.data.frame(cbind(tss, predict(basis, tss$time1))),
            names = c('id', 'time1', 'time2', paste0('basis', 1:deg))
  ))
}

# Using this StartStop and getting timevarying coxph...
TimeVarCox2 <- function(data, REs, basis, nK = 3, fixef.surv = c('cont', 'bin'),
                       survtime = 'survtime', status = 'status'){
  # Prepare data
  deg <- attr(basis, 'degree')
  ss <- ToStartStop2(data, basis)
  REs <- as.data.frame(REs); REs$id <- 1:nrow(REs)
  ss2 <- dplyr::left_join(ss, REs, 'id')
  ss2 <- dplyr::distinct(dplyr::left_join(ss2, data[, c('id', fixef.surv, survtime, status)], 'id'))
  
  # Create \gamma_k variables
  gammaK <- matrix(NA, nrow(ss2), nK)
  colnames(gammaK) <- paste0('gamma_', 1:nK)
  
  for(k in 1:nK){
    re.ii <- grepl(paste0('^b', k, '\\_'), x = colnames(ss2))
    re.k <- ss2[,re.ii]
    bases <- ss2[, grepl('^basis', colnames(ss2))]
    inter <- re.k[, grepl('\\(Intercept\\)' , colnames(re.k))]
    gammaK[, k] <- inter + rowSums(re.k[, !grepl('\\(Intercept\\)' , colnames(re.k))] * bases)
  }
  
  # And join on ...
  ss3 <- cbind(ss2, gammaK)
  # Update this to deal with ties too?
  ss3$status2 <- ifelse(ss3$survtime == ss3$time2, ss3$status, 0)
  
  # Time Varying coxph
  # Formula
  timevar.formula <- as.formula(
    paste0('Surv(time1, time2, status2) ~ ', paste0(fixef.surv, collapse = ' + '), ' + ', paste0('gamma_', 1:nK, collapse = ' + '))
  )
  ph <- coxph(timevar.formula, data = ss3)
  
  list(inits = coef(ph), l0.init = coxph.detail(ph)$haz, ph = ph)
}

