#' #######
#' Simulates data under a the following joint modelling specification:
#'   Longitudinal process: Where marginal distribution of Y given random effects is
#'                             assumed to belong to a member of the exponential family
#'   Survival process: Survival sub-model which is then linked with the above longitudinal process by some
#'                        shared random effects, which are assumed to be Gaussian.
#'                        
#'  Family: Can be "gaussian"; "binomial"; "poisson"; "negative.binomial".
#'  Must be a list of K families, and beta the corresponding fixed effects (intercept, time, standard normal, binary).
#'  The random effects structure is ALWAYS of the form intercept and slope.
#'  In cases where multiple of the same families are used, they are simulated with the same corresponding dispersion parameter.
#'     (maybe fix in future -- hard to think of immediate solution!)
#'  Supplied covariance matrix must be of dimension K*2, and gamma order K.
#'  Default is to generate K = 2 balanced longitudinal responses.
#' #######

simData <- function(n = 250, ntms = 10, fup = 5, family = list('gaussian', 'gaussian'), disp = NULL, 
                    beta = rbind(c(1, 0.10, 0.33, -0.50), c(1, 0.10, 0.33, -0.50)), D = NULL, var.e = 0.16,
                    gamma = c(0.5, -0.5), zeta = c(0.05, -0.30), theta = c(-4, 0.2), cens.rate = exp(-3.5),
                    random.formula = NULL){
  
  #' Checks --------------
  # Check family is valid option
  family <- lapply(family, function(f){
    if("function"%in%class(f)) f <- f()$family
    if(!f%in%c("gaussian", "binomial", "poisson", "negative.binomial")) stop('Family must be one of "gaussian", "binomial", "poisson" or "negative.binomial"')
    f
  })
  # Check dispersion supplied if neg. binom.
  family <- lapply(family, function(f) family <- match.arg(f, c('gaussian', 'binomial', 'poisson', 'negative.binomial'), several.ok = F))
  # Checks wrt the fixed effects and dispersion parameters
  funlist <- unlist(family)
  num.nb <- length(which(funlist == 'negative.binomial'))
  num.ga <- length(which(funlist == 'gaussian'))
  # Fixed effects
  K <- nrow(beta)
  if(K != length(funlist)) stop('Incorrect number of families provided and/or beta terms.')
  if(K != length(gamma)) stop('Incorrect number of association parameters provided wrt beta terms.')
  if(is.null(D)) D <- diag(K * 2) else D <- D
  if(!is.null(random.formula)){
    if(length(random.formula) != K) stop('Provided random formulas must be a list of length K.')
  }else{
    if((K*2) != ncol(D)) stop('Incorrect dimension on supplied covariance matrix D, do you need to supply random.formula?')
  }
  
  # Check covariance matrix D is positive semi-definite.
  if(any(eigen(D)$value < 0) || det(D) <= 0) stop('D must be positive semi-definite')
  
  #' Necessary parameters & data generation ----
  time <- seq(0, fup, length.out = ntms); tau <- fup + 0.1
  cont <- rnorm(n); bin <- rbinom(n, 1, 0.5)  # Continuous and binary covariates.
  
  df <- data.frame(id = rep(1:n, each = ntms),
                   time = rep(time, n),
                   cont = rep(cont, each = ntms),
                   bin = rep(bin, each = ntms))
  
  #' Design matrices, used across ALL responses ----
  X <- model.matrix(~ time + cont + bin, data = df)
  if(!is.null(random.formula)){
    Z <- lapply(random.formula, function(x) model.matrix(x, data = df))
  }else{
    Z <- replicate(K, model.matrix(~ 1 + time, data = df), simplify = F) # assume all intslopes.
  }

  #' Linear predictor & response generation ----
  if(!is.null(random.formula)){
    q <- ncol(do.call(cbind, lapply(Z, head)))
    b <- MASS::mvrnorm(n, mu = rep(0, q), Sigma = D)
    b.inds <- split(1:q, do.call(c, sapply(1:K, function(k) rep(k, ncol(Z[[k]])))))
  }else{
    b <- MASS::mvrnorm(n, mu = rep(0, K * 2), Sigma = D)
    b.inds <- split(1:(2*K), rep(1:K, each = 2)) 
  }
  Y <- sapply(1:K, function(k){
    f <- family[[k]]; Zk <- Z[[k]]
    betak <- beta[k,,drop=T]
    etak <- X %*% betak + rowSums(Zk * b[df$id, b.inds[[k]]])
    switch(f, 
           gaussian = Y <- rnorm(n * ntms, mean = etak, sd = sqrt(var.e)),
           binomial = Y <- rbinom(n * ntms, 1, plogis(etak)),
           poisson = Y <- rpois(n * ntms, exp(etak)),
           negative.binomial = Y <- MASS::rnegbin(n * ntms, mu = exp(etak), theta = disp)
    )
    Y
  })
  colnames(Y) <- paste0('Y.', 1:K)
  
  df <- cbind(df, Y)
  print(df)
  
  #' Survival ----
  theta0 <- theta[1]; theta1 <- theta[2]
  Keta <- cbind(cont, bin) %*% zeta
  U <- runif(n)
  if(!is.null(random.formula)){
    zz <- lapply(Z, colnames)
    ints <- which(grepl('\\(Intercept\\)', do.call(c, zz))); slopes <- which(grepl('time', do.call(c, zz)))
    if(any(!grepl('\\(Intercept\\)|time', do.call(c, zz)))) message('Warning: Only intercept-only or intercept-ands-slope random.formulas allowed.')
    # Which of the K responses have intercept/slope
    ints.gamma <- do.call(c, lapply(1:K, function(k){
      if(any(grepl('\\(Intercept\\)', colnames(Z[[k]])))) return(k) else return(NULL)
    }))
    slopes.gamma <- do.call(c, lapply(1:K, function(k){
      if(any(grepl('time', colnames(Z[[k]])))) return(k) else return(NULL)
    }))
  }else{
    ints <- seq(1,2*K,by=2); slopes <- seq(2, 2*K, by = 2) # intslope on all.
    ints.gamma <- slopes.gamma <- 1:K    # Each K has both intercept and slope.
  }
  
  b0 <- b[, ints, drop = F]; b1 <- b[, slopes, drop = F]
  print(b0); print(b1)
  print(gamma[ints.gamma]);print(gamma[slopes.gamma])
  # Generate survival times (Austin et al 2012)
  denom <- theta1 + b1 %*% gamma[slopes.gamma]
  rhs <- (theta1 + b1 %*% gamma[slopes.gamma]) * log(U)/(exp(theta0 + Keta + b0 %*% gamma[ints.gamma]))
  t <- suppressWarnings(log(1-rhs)/denom)
  t[is.nan(t)] <- tau
  
  # Collect survival times, and generate cenors times.
  cens.time <- rexp(n, cens.rate)
  survtime <- pmin(t, cens.time)
  survtime[survtime >= tau] <- tau
  
  # Status flag
  status <- rep(1, n)
  is.censored <- cens.time < survtime
  status[which(survtime == tau | is.censored | survtime == cens.time)] <- 0 # Failure flag
  
  #' Output Dataframes ----
  surv.data <- data.frame(id = 1:n, survtime, status)
  long.data <- df
  
  out.data <- dplyr::left_join(df, surv.data, by = 'id')
  out.data <- out.data[out.data$time < out.data$survtime, ]
  message(round(100 * sum(surv.data$status)/n), '% failure rate')
  
  list(data =  out.data, 
       surv.data =  dplyr::distinct(out.data, id, survtime, status, cont, bin))
}
