#' #####
#' Dynamic survival predictions
#' #####

dynSurv <- function(fit, data, id, u = NULL, b.method = 'normal', nsim = 200){
  # Checks
  if(!b.method %in% c('normal', 'MH')) stop('\nb.method should be either normal or MH.\n')
  tmax <- max(unique(data[data$status == 1, 'survtime']))
  if(any(u > tmax)) stop("Can't extrapolate beyond last failure time.")
  newdata <- data[data$id == id, ]
  if(is.null(u)){
    u <- newdata$time + 1 # all but time  = 0 if not specified
    if(any(u > tmax)) u <- u[!which(u > tmax)] # and ensure those after T_{max} aren't included
  }
  
  pi <- setNames(vector('list', length = length(u)), paste0('u = ', u))
  for(uu in 1:length(u)){
    newdata2 <- newdata[newdata$time < u[uu], ]
    pt <- prepdata(newdata2, id = id, fit = fit)
    pu <- prepdata(newdata2, id = id, u = u[uu], fit = fit)
    
    if(b.method == 'MH'){
      b <- pt$b; Sigmai.prop <- pt$S #* 2
    }
    
    pi.store <- numeric(nsim)
    pb <- utils::txtProgressBar(max = nsim, style = 3)
    for(i in 1:nsim){
      O <- Omega.draw(fit)
      if(b.method == 'MH'){
        mh <- b.mh(O, Sigmai.prop, pt$b, pt)
        b <- mh$b
      }else{
        b <- b.draw(pt$b, 
                    pt$long$Xt, pt$long$Yt, pt$long$Zt,
                    O$beta, O$D,
                    pt$surv$Delta, pt$surv$K, pt$surv$Fi, pt$surv$l0i, pt$surv$KK.t,
                    pt$surv$Fu.t, pt$surv$l0u.t, O$gamma, O$eta)
        b <- b$b
      }
      
      pi.store[i] <- S(b, O, pu$surv) / S(b, O, pt$surv)
      utils::setTxtProgressBar(pb, i)
      
    }
    pi[[uu]] <- pi.store
  }
  cat('\n\n')
  return(lapply(pi, quantile, probs = c(.500, 0.025, .975)))
  
}