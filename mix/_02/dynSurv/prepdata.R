#' #######
#' Functions to prepare data for dynamic predictions
#' #######


# Prepare survival objects ------------------------------------------------
prep.surv <- function(data, l0, u = NULL, smoothed = F){
  # data: i.dat-type object (i.e. just one subject)
  # l0: final baseline hazard from fitted model (matrix, second column is hazard, first is fts)
  # u: if not null, then outputs both objects for u and t <= u for the subject.
  #    if null (default) then this is same as surv.mod with reduced output.

  udat <- dplyr::distinct(data[,c('cont', 'bin', 'survtime', 'status')])
  Delta <- unique(data$status)
  survtime <- unique(data$survtime)
  ft <- l0[,1]
  # print(c(Delta, survtime))
  if(smoothed) l0 <- predict(lm(l0 ~ splines::bs(ft, df = 4)))
  
  # Failure-time design objects
  Fi <- do.call(c, replicate(3, c(1, survtime), simplify = F))
  Fu.t <- cbind(1, l0[ft <= max(data$time), 1])
  if(length(which(ft <= max(data$time))) == 0) Fu.t <- cbind(1, 0)
  if(!is.null(u)){
    Fu.t <- cbind(1, l0[ft <= u, 1])
    # if(length(which(ft <= u)) == 0) Fu.t <- cbind(1, 0)
    # Fu.u <- cbind(1, l0[ft <= u, 1])
  }
  #print(survtime) 
  l0u.t <- l0[ft <= max(data$time), 2]
  if(length(which(ft <= max(data$time))) == 0) l0u.t <- 0
  if(!is.null(u)){
    l0u.t <- l0[ft <= u, 2]
    # l0u.u <- l0[ft <= u, 2]
  }
  if(Delta == 0) l0i <- 0 else l0i <- l0[which(ft == survtime), 2]
  # K, design matrices
  K <- unname(cbind(udat$cont, udat$bin))
  KK.t <- apply(K, 2, rep, nrow(Fu.t))
  if('numeric' %in% class(KK.t)) KK.t <- t(as.matrix(KK.t))
  # KK.u <- apply(K, 2, rep, nrow(Fu.u))
  # if('numeric' %in% class(KK.u)) KK.u <- t(as.matrix(KK.u))
  
  # output
  return(list(
    u = u,
    Delta = Delta, # indicator
    K = K, KK.t = KK.t, #KK.u = KK.u,  # K
    Fu.t = Fu.t, Fi = Fi, # F_(design)
    l0u.t = l0u.t, l0i = l0i # Hazard objects //
  ))
}


# Prepare longitudinal data -----------------------------------------------
prep.long <- function(data, u = NULL){
  # data: i.dat-type object (i.e. just one subject)
  # fit: Fitted model object containing coeffs
  # u: if not null, then outputs both objects for u and t <= u for the subject.
  #    if null (default) then gives full data matrices
  if(is.null(u)){
    Xt <- model.matrix(~ time + cont + bin, data)
    Yt <- cbind(data$Y.1, data$Y.2, data$Y.3)
    Zt <- model.matrix(~ time, data)
  }else{
    datat <- data[data$time <= u, ] # this a very bad way of doing this...
    Xt <- model.matrix(~ time + cont + bin, datat)
    Yt <- cbind(datat$Y.1, datat$Y.2, datat$Y.3)
    Zt <- model.matrix(~ time, datat)
    # datau <- data[data$time <= u, ] 
    # Xu <- model.matrix(~ time + cont + bin, datau)
    # Yu <- datau$Y
    # Zu <- model.matrix(~ time, datau)
  }
  return(list(
    Xt = Xt, Yt = Yt, Zt = Zt#,
    # Xu = Xu, Yu = Yu, Zu = Zu
  ))
  
}


# Prepare data ------------------------------------------------------------
prepdata <- function(data, id, u = NULL, fit, smoothed = F){
  # data = dataset containing ALL items
  # id = id of subject we want to get dynamic predictions for
  # u = candidate time we want dynpreds for
  
  i.dat <- data#[data$id == id, ]
  
  long <- prep.long(i.dat, u)
  surv <- prep.surv(i.dat, fit$hazard, u, smoothed)
  
  b <- fit$RE[id, ] # think of a better way to do this !
  
  # Establish if NB fit or not.
  if(length(fit$coeffs$theta) == 0){
    nb <- F; theta <- 0
  }else{
    nb <- T; theta <- fit$coeffs$theta
  }
  
  if(is.null(u)){
    
    V <- diag(fit$coeffs$var.e, nrow = nrow(long$Y), ncol = nrow(long$Y)) # Create matrix of residual variance
    S <- solve(joint_density_sdb(
      b, long$Xt, long$Zt, fit$coeffs$beta, V, fit$coeffs$D,
      long$Yt[, 1], long$Yt[, 2], long$Yt[, 3],
      nb, theta, 
      surv$Delta, surv$K, surv$Fi, surv$l0i, surv$KK.t, surv$Fu.t,
      surv$l0u.t, rep(fit$coeffs$gamma, each = 2), fit$coeffs$eta, 1e-3
    ))
  }
  
  
  list(long = long, surv = surv, b = b, S = S)
  
}
