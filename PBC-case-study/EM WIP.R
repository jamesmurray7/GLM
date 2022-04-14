vech <- function(x) x[lower.tri(x, T)]


EMupdate <- function(log.likelihood, score.eta, vare.update = NULL, theta.update = NULL,
                     mu, Omega, X, Y, Z, b, Fi, Fu, l0i, l0u, Delta, l0, sv, surv, w, v){
  
}

EM <- function(long.formula, surv.formula, data, family, control = list()){
  
  formulas <- parseFormula(long.formula)
  surv <- parseCoxph(surv.formula, data)
  n <- nrow(surv$survdata)
  
  
  #' Initial conditons ----
  if(!is.null(control$dispformula)) dispformula <- control$dispformula else dispformula <- NULL
  inits.long <- Longit.inits(long.formula, data, family, dispformula = dispformula)
  inits.surv <- TimeVarCox(data, inits.long$b, surv$ph)
  
  # Longitudinal parameters
  beta <- inits.long$beta.init
  D <- inits.long$D.init
  theta <- inits.long$theta.init
  var.e <- inits.long$var.e.init
  b <- lapply(1:n, function(i) inits.long$b[i, ])
  # Survival parameters
  zeta <- inits.surv$inits[match(colnames(surv$ph$x), names(inits.surv$inits))]
  names(zeta) <- paste0('zeta_', names(zeta))
  gamma <- inits.surv$inits[grepl('gamma', names(inits.surv$inits))]
  
  #' Longitudinal and survival data objects ----
  sv <- surv.mod(surv$ph, surv$survdata, formulas, inits.surv$l0.init)
  dmats <- createDataMatrices(data, formulas)
  X <- dmats$X; Y <- dmats$Y; Z <- dmats$Z # Longitudinal data matrices
  Fi <- sv$Fi; Fu <- sv$Fu; l0i <- sv$l0i; l0u <- sv$l0u; Delta <- surv$Delta # survival
  l0 <- sv$l0
  
  
  #' Parameter vector ----
  Omega <- list(D=D, beta = beta, vare = var.e, theta = theta, zeta = zeta, gamma = gamma)
  params <- c(setNames(vech(D), paste0('D[', apply(which(lower.tri(D, T), arr.ind = T), 1, paste, collapse = ','), ']')),
              beta, vare = var.e, theta = theta, zeta, gamma)
  
  
  #' Gauss-Hermite Quadrature ----
  if(!is.null(control$gh.nodes)) gh <- control$gh.nodes else gh <- 3
  if(!is.null(control$gh.sigma)) .sigma <- control$gh.sigma else .sigma <- 1
  
  GH <- statmod::gauss.quad.prob(gh, 'normal', sigma = .sigma)
  w <- GH$w; v <- GH$n
  
  #' Assign joint density (& derivs) along with S/Hbeta (and updates to theta and vare if appropriate) -----------
  if("function"%in%class(family)) family <- family()$family # step to ensure non-quoted arguments don't throw error.
  family <- match.arg(family, c('gaussian', 'binomial', 'poisson', 'negative.binomial'), several.ok = F)
  
  switch(family,  ## THIS WIP !!!!!!!!!!!!!!!!!!!!!
         gaussian = {
           log.likelihood <- gauss_ll
           score.eta <- Score_eta_gauss
           vare.update <- vare_update
           mu <- identity
         },
         binomial = {
           log.likelihood <- binomial_ll
           score.eta <- Score_eta_binomial
           mu <- plogis
         },
         poisson = {
           log.likelihood <- poisson_ll
           score.eta <- Score_eta_poisson
           mu <- exp
         },
         negative.binomial = {
           log.likelihood <- negbin_ll
           score.eta <- Score_eta_negbin
           mu <- exp
           theta.update <- theta_update
         }
  )
  
  
  
}