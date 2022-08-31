ff2 <- function(delta, Y, G, mu, summax){
  nu <- exp(G %*% delta)
  VY <- expected_variance(mu, nu, summax)
  sum(c(var(Y)) - VY)
}

# Determine interval in the 1-d case.
determine.interval <- function(mu, summax){
  m <- length(mu)
  lower <- -3; upper <- 3 # Nominally choose lower and upper values of interval to check over.
  test.lower <- lambda_appx(mu, rep(exp(lower), m), summax)
  test.upper <- lambda_appx(mu, rep(exp(upper), m), summax)
  
  flag.lower <- any(is.nan(test.lower)); flag.upper <- any(is.nan(test.upper))
  while(flag.lower){
    lower <- lower + 0.2
    flag.lower <- any(is.nan(lambda_appx(mu, rep(exp(lower), m), summax)))
  }
  while(flag.upper){
    upper <- upper - 0.2
    flag.upper <- any(is.nan(lambda_appx(mu, rep(exp(upper), m), summax)))
  }
  c(lower = lower, upper = upper)
}


delta.uniroot <- function(Y, G, mu, summax){
  .interval <- determine.interval(mu, summax)
  tryCatch(uniroot(ff2, interval = c(.interval[1], .interval[2]), Y = Y, G = G, mu = mu, summax = summax, extendInt = 'yes')$root,
           error = function(e) NA)
}

find.deltas <- function(Ylist, Glist, mulist, summax, verbose = F, min.profile.length = 1){
  candidateY <- sapply(Ylist, length) > min.profile.length
  numY <- sum(candidateY); inds <- which(candidateY)
  out <- numeric(numY); p <- 1
  for(i in inds){
    out[p] <- delta.uniroot(Ylist[[i]], Glist[[i]], mulist[[i]], summax)
    if(verbose) cat(sprintf('%d/%d\r', p, numY))
    p <- p + 1
  }
  out
}

# Optim version ----
ff3 <- function(delta, Y, G, mu, summax){
  nu <- exp(G %*% delta)
  VY <- expected_variance(mu, nu, summax)
  abs(sum(c(var(Y)) - VY))
}

# Determine interval for the n-dimensional case...
# ... (wip)
determine.interval.nd <- function(Y, G, mu, summax){
  g <- ncol(G)
  
  # we want g permutations of seq(-3,3,.2)
  candidates <- seq(-3, 3, .2)
  candidates.mat <- structure(t(combn(candidates, g)),
                              dimnames = list(as.character(1:choose(length(candidates),g)),
                                              paste0('delta', 1:g))
  )
  
  vars <- sort(apply(candidates.mat, 1, function(i){
    ff2(c(i), Y[[1]], G[[1]], mus[[1]], summax)
  }))
}

delta.optim <- function(Y, G, mu, summax){
  g <- ncol(G)
  if(g == 1){
    # .interval <- determine.interval(mu, summax)
    out <- tryCatch(optim(c(0), ff3, NULL,
                          Y = Y, G = G, mu = mu, summax = summax,
                          method = 'Brent', lower = -2, upper = 2)$par,
                    error = function(e) NA)
  }else{
    out <- tryCatch(optim(rep(0, g), ff3, NULL,
                          Y = Y, G = G, mu = mu, summax = summax,
                          method = 'L-BFGS-B', lower = rep(-2, g), upper = rep(2, g))$par,
                    error = function(e) NA)
  }
  out
}

find.deltas.optim <- function(Ylist, Glist, mulist, summax, verbose = F, min.profile.length = 1){
  candidateY <- sapply(Ylist, length) > min.profile.length
  numY <- sum(candidateY); inds <- which(candidateY)
  out <- numeric(numY); p <- 1
  for(i in inds){
    out[p] <- delta.optim(Ylist[[i]], Glist[[i]], mulist[[i]], summax)
    if(verbose) cat(sprintf('%d/%d\r', p, numY))
    p <- p + 1
  }
  out
}

# mapply(function(Y,G,mu){
#   g <- ncol(G)
#   x <- tryCatch(optim(rep(0, g), ff3, NULL,
#              Y = Y, G = G, mu = mu, summax = summax,
#              method = 'L-BFGS-B', lower = rep(-2, g), upper = rep(2, g)),
#              error = function(e) NULL)
#   if(!is.null(x)) return(list(x$par, x$value, x$conver)) else return(rep(NA, 3))
# }, Y = Y, G = G, mu = mus, SIMPLIFY = F) -> fits
# 
# obj <- do.call(c, lapply(fits, el, 2))
# pars <- do.call(rbind, lapply(fits, el, 1))
# apply(pars[obj < 1e-3 & !is.na(obj),],2,median)
