setwd('~/Documents/GLMM/zip')
sourceCpp('zip.cpp')
source('fresh/_Functions.R')

#' ----
#' Set out true parameters
#' ----
beta <- c(1.5, 0.05, 0.33, 0.50)
alpha <- c(-0.5, 0.25)
D <- diag(c(.5^2, .15^2))
n <- 250
ntms <- 15
data <- replicate(100, 
                  simData_zip(n, ntms, beta, alpha, D), simplify = F)

#' ----
#' Define EM step
#' ----
EMupdateCpp <- function(b, Y, X, Z, Xz, Zz, 
                        beta, D, alpha, gh, indzi = 2){
  
  b.hat <- mapply(function(b, Y, X, Z, Xz, Zz){
    ucminf::ucminf(b, b_logdensity, b_score, Y, X, Z, Xz, Zz, 
                   beta, alpha, D, indzi)$par
  }, b = b, Y = Y, X = X, Z = Z, Xz = Xzi, Zz = Zzi, SIMPLIFY = F)
  
  Sigmai <- mapply(function(b, Y, X, Z, Xz, Zz){
    solve(fd_b(b, Y, X, Z, Xz, Zz, beta, alpha, D, indzi, 1e-4))
  }, b = b.hat, Y = Y, X = X, Z = Z, Xz = Xzi, Zz = Zzi, SIMPLIFY = F)
  
  Drhs <- mapply(function(b, S){
    S + tcrossprod(b)
  }, b = b.hat, S = Sigmai, SIMPLIFY = F)
  
  ba <- mapply(function(b, Y, X, Z, Xz, Zz, Sigmai){
    beta_alpha_new(b, Y, X, Z, Xz, Zz, beta, D, alpha, Sigmai, gh)
  }, b = b.hat, Y = Y, X = X, Z = Z, Xz = Xzi, Zz = Zzi, Sigmai = Sigmai, SIMPLIFY = F)
  
  # Updates
  Sbeta <- Reduce('+', lapply(ba, '[[', 1))
  Hbeta <- Reduce('+', lapply(ba, '[[', 3))
  Salpha <- Reduce('+', lapply(ba, '[[', 2))
  Halpha <- Reduce('+', lapply(ba, '[[', 4))
  
  D.new <- Reduce('+', Drhs)/n
  beta.new <- beta - solve(Hbeta, Sbeta)
  alpha.new <- alpha - solve(Halpha, Salpha)
  
  return(list(
    D.new = D.new, beta.new = beta.new, alpha.new = alpha.new, b.hat = b.hat
  ))
}

#' ----
#' Define mini function
#' ----
vech <- function(x) x[lower.tri(x,diag=T)]
EM <- function(data, indzi = 2, gh){
  start <- proc.time()[3]
  inits.long <- Longit.inits(data)
  beta <- inits.long$beta.init
  alpha <- inits.long$alpha.init
  D <- inits.long$D.init
  b <- Ranefs(inits.long); n <- nrow(b)
  b <- lapply(1:n, function(i) b[i, ])
  # Get data matrices...
  X <- Y <- Z <- Xz <- Zz <- list()
  for(i in 1:n){
    i.dat <- data[data$id == i, ]
    X[[i]] <- model.matrix(~time + cont + bin, i.dat)
    Xz[[i]] <- model.matrix(~time, i.dat)
    Z[[i]] <- Zz[[i]] <- model.matrix(~1, i.dat)
    Y[[i]] <- i.dat$Y
  }          
  EMstart <- proc.time()[3]
  params <- c(vech(D), beta, alpha)
  diff <- 100; iter <- 0
  while(diff > tol){
    update <- EMupdateCpp(b, Y, X, Z, Xz, Zz, beta, D, alpha, gh = gh)
    params.new <- c(vech(update$D.new), update$beta.new, update$alpha)
    diff <- max(
      abs(params.new-params)/(abs(params)+1e-3)
    )
    message('\nIteration ', iter + 1)
    message('Maximum relative difference: ', round(diff, 5))
    b <- update$b.hat
    D <- update$D.new
    beta <- update$beta.new
    alpha <- update$alpha.new
    params <- c(vech(D), beta, alpha)
    iter <- iter + 1
  }
  EMend <- proc.time()[3]
  coeffs <- list(D = D, beta = beta, alpha = alpha, b = b)
  out <- list(coeffs = coeffs,
              inits = inits.long,
              EMtime = EMend-EMstart,
              iter = iter,
              totaltime = proc.time()[3] - start)
  out
}

#' ----
#' LOOP
#' ----
pb <- utils::txtProgressBar(max = 100, style = 3)
fits <- list()
for(m in 1:100){
  fits[[m]] <- tryCatch(suppressMessages(EM(data[[m]], 2, 3)), 
                        error = function(e) NULL)
  utils::setTxtProgressBar(pb, m)
}

