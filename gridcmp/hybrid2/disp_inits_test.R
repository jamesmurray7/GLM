#' Testing dispersion parameter initial estimates
#' on different true dispersion(s).

# Simulate some data ------------------------------------------------------
sets.8 <- replicate(100,
                    simData_joint2(n = 250, delta = c(.8, 0), 
                                   ntms = 10, theta = c(-2, .1), fup = 3,
                                   beta = c(2, -0.1, 0.1, -0.2), gamma = 0.6, zeta= c(0.0, -0.2),
                                   D = matrix(c(0.25, 0, 0, 0.00), 2, 2))$data,
                    simplify = F)

sets.3 <- replicate(100,
                    simData_joint2(n = 250, delta = c(.3, 0), 
                                   ntms = 10, theta = c(-2, .1), fup = 3,
                                   beta = c(2, -0.1, 0.1, -0.2), gamma = 0.6, zeta= c(0.0, -0.2),
                                   D = matrix(c(0.25, 0, 0, 0.00), 2, 2))$data,
                    simplify = F)

sets.minus5 <- replicate(100,
                    simData_joint2(n = 250, delta = c(-.5, 0), 
                                   ntms = 10, theta = c(-2, .1), fup = 3,
                                   beta = c(2, -0.1, 0.1, -0.2), gamma = 0.6, zeta= c(0.0, -0.2),
                                   D = matrix(c(0.25, 0, 0, 0.00), 2, 2))$data,
                    simplify = F)


long.formula <- Y ~ time + cont + bin + (1|id)
surv.formula <- Surv(survtime, status) ~ bin
disp.formula <- ~1

# Save data
save(sets.8, file = '/data/c0061461/disp-data-8.RData')
save(sets.3, file = '/data/c0061461/disp-data-3.RData')
save(sets.minus5, file = '/data/c0061461/disp-data-minus5.RData')

get.delta.inits <- function(data){
  #' Parsing formula objects ----
  formulas <- parseFormula(long.formula)
  N <- 250
  
  #' Initial conditions ----
  inits.long <- Longit.inits(long.formula, disp.formula, data)
  
  # Longitudinal parameters
  beta <- inits.long$beta.init
  D <- inits.long$D.init
  b <- lapply(1:N, function(i) inits.long$b[i, ])
  delta <- inits.long$delta.init
  
  #' Data objects ----
  dmats <- createDataMatrices(data, formulas, disp.formula)
  X <- dmats$X; Y <- dmats$Y; Z <- dmats$Z # Longitudinal data matrices
  lY <- lapply(Y, lfactorial)
  G <- dmats$G                             # Dispersion data matrix
  m <- sapply(Y, length)
  
  summax <- ceiling(max(sapply(Y, max)) * 2)
  mus <- mapply(function(X, Z, b) exp(X %*% beta + Z %*% b), X = X, Z = Z, b = b)
  
  a <- proc.time()[3]
  # raw <- find.deltas(Y, G, mus, summax, F)
  raw <- find.deltas.optim(Y, G, mus, summax, T)
  b <- proc.time()[3]
  list(
    subject.estimates = raw,
    median.estimate = median(raw, na.rm = T),
    IQR.estimates = IQR(raw, na.rm = T),
    sd.estimates = sd(raw, na.rm = T),
    time = round(b-a,3)
  )
  
}

load('/data/c0061461/disp-data-3.RData')
load('/data/c0061461/disp-data-8.RData')
load('/data/c0061461/disp-data-minus5.RData')

disp.inits.8 <- lapply(sets.8, get.delta.inits)
save(disp.inits.8, file = '/data/c0061461/disp-fits-8.RData')
disp.inits.3 <- lapply(sets.3, get.delta.inits)
save(disp.inits.3, file = '/data/c0061461/disp-fits-3.RData')
disp.inits.minus5 <- lapply(sets.minus5, get.delta.inits)
save(disp.inits.minus5, file = '/data/c0061461/disp-fits-minus5.RData')
