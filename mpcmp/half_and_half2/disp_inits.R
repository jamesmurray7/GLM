# Optim version ----
ff3 <- function(delta, Y, mu, summax){
  nu <- rep(exp(delta), length(Y))
  VY <- expected_variance(mu, nu, summax)
  abs(sum(c(var(Y)) - VY))
}

make.interval <- function(Y){
  disp <- as.double(var(Y)/mean(Y));
  if(disp > 1) # OD
    return(c(-10, 0))
  else
    return(c(0, 10))
}

delta.optim <- function(Y, mu, summax){
  bounds <- make.interval(Y)
  # out <- tryCatch(optim(c(mean(bounds)/5), ff3, NULL,
  #                       Y = Y, mu = mu, summax = summax,
  #                       method = 'Brent', lower = bounds[1], upper = bounds[2])$par,
  #                 error = function(e) NA)
  # out
  optim(c(mean(bounds)/5), ff3, NULL,
        Y = Y, mu = mu, summax = summax,
        method = 'Brent', lower = bounds[1], upper = bounds[2])$par
}


find.deltas.optim <- function(Ylist, mulist, summax, verbose = F, min.profile.length = 1){
  out <- numeric(length(Ylist))
  candidateY <- sapply(Ylist, function(y) length(unique(y))) > min.profile.length
  numY <- sum(candidateY); inds <- which(candidateY); p <- 1
  for(i in inds){
    out[i] <- this <- delta.optim(Ylist[[i]], mulist[[i]], summax[[i]])
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
  out <- numeric(length(Ylist))
  candidateY <- sapply(Ylist, function(y) length(unique(y))) > min.profile.length
  numY <- sum(candidateY); inds <- which(candidateY); p <- 1
  for(i in inds){
    out[i] <- delta.bobyqa(Ylist[[i]], mulist[[i]], summax[[i]])
    if(verbose) cat(sprintf('%d/%d\r', p, numY))
    p <- p + 1
  }
  out
}
