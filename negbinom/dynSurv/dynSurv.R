#' #####
#' Dynamic survival predictions
#' #####

dynSurv <- function(fit, data, id, u = NULL, nsim = 200, progress = T){
  # Checks
  tmax <- max(unique(data[data$status == 1, 'survtime']))
  if(any(u > tmax)) stop("Can't extrapolate beyond last failure time.")
  newdata <- data[data$id == id, ]
  if(is.null(u)){
    u <- newdata$time + 1 # all but time  = 0 if not specified
    if(any(u > tmax)) u <- u[!which(u > tmax)] # and ensure those after T_{max} aren't included
  }
  
  newdata2 <- newdata[newdata$time <= u[1], ]
  pt <- prepdata(newdata2, id = id, fit = fit)
  
  u <- u[-1] # Don't want to find dynsurv from u = t (Tstart)
  
  pi <- setNames(vector('list', length = length(u)), paste0('u = ', u))
  for(uu in 1:length(u)){
    pu <- prepdata(newdata2, id = id, u = u[uu], fit = fit)

    pi.store <- numeric(nsim)
    if(progress) pb <- utils::txtProgressBar(max = nsim, style = 3)
    for(i in 1:nsim){
      
      error.flag <- T
      while(error.flag){
        error.flag.pi <- T
        O <- Omega.draw(fit)
        b <- tryCatch(b.draw(pt$b, 
                             pt$long$Xt, pt$long$Yt, pt$long$Zt,
                             O$beta, O$theta, O$D,
                             pt$surv$Delta, pt$surv$K, pt$surv$Fi, pt$surv$l0i, pt$surv$KK.t,
                             pt$surv$Fu.t, pt$surv$l0u.t, O$gamma, O$eta),
                      error = function(e) NA)
        error.flag.b <- any(is.na(b))
        if(!error.flag.b){
          b <- b$b
          pi.store[i] <- S_Cpp2(pt$surv, pu$surv, rep(O$gamma, each = 2), O$eta, b)
          error.flag.pi <- is.nan(pi.store[i])
        }
        error.flag <- any(error.flag.b, error.flag.pi)
      }

      if(progress) utils::setTxtProgressBar(pb, i)
      
    }
    pi[[uu]] <- pi.store
  }
  # cat('\n\n')
  return(lapply(pi, quantile, probs = c(0.500, 0.025, 0.975), na.rm = T))
}

# ROC and AUC -------------------------------------------------------------
ROC <- function(fit, data, Tstart, Delta, what = 'lowest', ...){
  # fit: Fitted object
  # data: Full data
  # Tstart: Time from which all survival probabilities are to be calculated from
  # Delta: Tie interval to check for failures in.
  # what: 'lowest': Lowest p determines threshold comparison; 'last' is last time's p in window.
  # ...: Additional arguments to dynSurv()
  if(!what%in%c('last', 'lowest')) stop("'What' should be either 'last' or 'lowest'.")
  tmax <- max(data[data$status == 1, 'survtime'])
  # Set-out new data
  newdata <- data[data$survtime > Tstart, ] # people ALIVE AND IN STUDY at Tstart
  alive.ids <- unique(newdata$id)
  n.alive <- length(alive.ids)
  # Failure times to sample from and sampling window
  fts <- sort(unique(newdata[newdata$status == 1, 'survtime', drop = T]))
  window <- c(Tstart + 1e-6, Tstart + 1e-6 + Delta)
  if(window[2] > tmax) window[2] <- tmax
  # Candidate failure times
  candidate.u <- c(Tstart, fts[fts >= window[1] & fts <= window[2]])
  probs <- setNames(vector('list', length = length(alive.ids)), 
                    paste0('id ', alive.ids))
  # Loop over ids and failure times
  pb <- utils::txtProgressBar(max = length(alive.ids), style = 3)
  for(i in seq_along(alive.ids)){
    ds <- dynSurv(fit, newdata, id = alive.ids[i], u = candidate.u, progress = F,  ...)
    probs[[i]] <- do.call(rbind, ds)#[-1, ] 
    utils::setTxtProgressBar(pb, i)
  }
  close(pb)
  
  # Working out whether individuals failed in the window
  events <- with(newdata, tapply(survtime, id, function(x){
    x >= window[1] && x <= window[2] # Any event in the window
  }))
  status <- as.logical(with(newdata, tapply(status, id, unique))) # Check if they failed
  event <- status & events # Check if they failed in window.
  
  # Obtaining conditional probabilities for sample alive at Tstart
  infodf <- lapply(alive.ids, function(x){
    p <- as.data.frame(probs[[paste0('id ', x)]])
    p$id <- x
    p 
  })
  # pi(u|t)
  if(what == 'last') pi <- with(do.call(rbind, infodf), tapply(`50%`, id, tail, 1))
  else if(what == 'lowest') pi <- with(do.call(rbind, infodf), tapply(`50%`, id, min))
  
  # Defining threshold and deriving metrics
  t <- seq(0, 1, length = 101)
  simfail <- structure(outer(pi, t, '<'),
                       dimnames = list(names(pi), paste0('t: ', t)))
  
  TP <- colSums(c(event) * simfail)   # True positives
  FN <- sum(event) - TP               # False negatives
  TPR <- TP/(TP + FN)                 # True positive rate (sensitivity)
  FP <- colSums(c(!event) * simfail)  # False positives
  TN <- sum(!event) - FP              # True negatives
  FPR <- FP / (FP + TN)               # False positive rate (1-specificity)
  
  # sanity checks
  if(!identical(TPR, TP/sum(event))) stop('Something wrong: TP + FN != sum(event)')
  if(!identical(FP / (FP + TN)  , FP/(sum(!event)))) stop('Something wrong: FP + TN != sum(!event)')
  
  # Making nice data.frame
  out.df <- data.frame(threshold = t,
                       TP = TP, TN = TN, FP = FP, FN = FN,
                       TPR = TPR, FPR = FPR)
  
  # Removing duplicated rows
  flip.out.df <- out.df[match(out.df$threshold, sort(out.df$threshold, decreasing = T)), ]
  metrics.full <- out.df
  metrics <- flip.out.df[rowSums(apply(flip.out.df[, c('TPR', 'FPR')], 2, duplicated))!=2,]
  
  # Deriving cutoff
  cutoff <- with(metrics, which.max(TPR-FPR))
  
  return(list(
    metrics = metrics,
    metrics.full = metrics.full,
    num.events = sum(events),
    num.ids = alive.ids,
    Tstart = Tstart, Delta = Delta,
    what = what, cutoff = cutoff
  ))
}

# Function to plot this object
plotROC <- function(ROC, cutoff = F, legend = F){
  TPR <- ROC$metrics$TPR; FPR <- ROC$metrics$FPR;
  plot(FPR, TPR,
       xlab = '1 - Specficity', ylab = 'Sensitivity',
       main = paste0('ROC curve for interval (', ROC$Tstart, ', ', ROC$Tstart + ROC$Delta, ']'),
       type = 'l')
  abline(0, 1, lty = 3)
  if(cutoff) abline(v = ROC$metrics[ROC$cutoff, 'FPR'], lty = 3, col = 'red')
  if(legend){
    legend('bottomright', 
           paste0(length(ROC$num.ids), ' at risk; ', ROC$num.events, ' failures in interval.\n',
                  'AUC: ', round(AUC(ROC), 3)),
           bty = 'n', cex = .75)
  }
  invisible()
}

# Determining area under the ROC curve.
AUC <- function(ROC){
  TPR <- rev(ROC$metrics$TPR); FPR <- rev(ROC$metrics$FPR);
  auc <- sum(0.5 * diff(FPR) * (TPR[-1] + TPR[-length(TPR)]), na.rm = TRUE)
  auc
}

