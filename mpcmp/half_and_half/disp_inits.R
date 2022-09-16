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
ff3 <- function(delta, Y, mu, summax){
  nu <- rep(exp(delta), length(Y))
  VY <- expected_variance(mu, nu, summax)
  abs(sum(c(var(Y)) - VY))
}

delta.optim <- function(Y, mu, summax){
  out <- tryCatch(optim(c(0), ff3, NULL,
                        Y = Y, mu = mu, summax = summax,
                        method = 'Brent', lower = -2, upper = 2)$par,
                  error = function(e) NA)
  out
}

find.deltas.optim <- function(Ylist, mulist, summax, verbose = F, min.profile.length = 1){
  out <- numeric(length(Y))
  candidateY <- sapply(Ylist, length) > min.profile.length
  numY <- sum(candidateY); inds <- which(candidateY)
  out <- numeric(numY); p <- 1
  for(i in inds){
    out[i] <- delta.optim(Ylist[[i]], mulist[[i]], summax)
    if(verbose) cat(sprintf('%d/%d\r', p, numY))
    p <- p + 1
  }
  out
}


# Bounded non-linear optimisation by quadratic approximation of the objective function.
delta.bobyqa <- function(Y, mu, summax){
  out <- tryCatch(nloptr::bobyqa(c(0), ff3,
                                 Y = Y, mu = mu, summax = summax)$par,
                  error = function(e) NA)
  out
}

find.deltas.bobyqa <- function(Ylist, mulist, summax, verbose = F, min.profile.length = 1){
  out <- numeric(length(Y))
  candidateY <- sapply(Ylist, function(y) length(unique(y))) > min.profile.length
  numY <- sum(candidateY); inds <- which(candidateY); p <- 1
  for(i in inds){
    out[i] <- delta.bobyqa(Ylist[[i]], mulist[[i]], summax)
    if(verbose) cat(sprintf('%d/%d\r', p, numY))
    p <- p + 1
  }
  out
}
