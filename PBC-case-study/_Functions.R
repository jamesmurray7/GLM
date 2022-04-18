#' ######
#' _Functions.R
#'   Collection of functions to help parse formula inputs
#' ######

# Parsing input formula ---------------------------------------------------
# (for longitudinal part), initial conditions found in inits.R
parseFormula <- function(formula){
  if(class(formula) == "formula") f <- as.character(formula) else stop('Argument "formula" not of class "formula".')
  
  response <- f[2]; rhs <- f[3]
  
  #' Split into fixed and random -----
  rhs <- gsub('\\s', '', rhs)
  split <- el(strsplit(rhs, '\\+\\(')) 
  # Fixed 
  fixed <- as.formula(paste0(' ~ ', split[1]))
  # Random
  full.random <- gsub('\\)', '', split[2])
  random.split <- el(strsplit(full.random, '\\|'))
  random <- as.formula(paste0('~', gsub('1\\+', '', random.split[1])))
  
  list(
    response = response, 
    fixed = fixed,
    random = random,
    random.by = random.split[2] # dont think this will ever get used.
  )
}


# Longitudinal objects ----------------------------------------------------
createDataMatrices <- function(data, formulas){
  fixed <- formulas$fixed
  response <- formulas$response
  random <- formulas$random
  
  n <- length(unique(data$id))
  X <- Z <- Y <- vector('list', n)
  for(i in 1:n){
    i.dat <- subset(data, id == i)
    X[[i]] <- model.matrix(fixed, i.dat)
    Y[[i]] <- matrix(i.dat[, response], nc = 1)
    Z[[i]] <- model.matrix(random, i.dat)
  }
  
  list(
    X = X, Y = Y, Z = Z
  )
}

# Survival objects --------------------------------------------------------
# Parsing the coxph object/formula (formula as of 14/4/22) and constructing 
# all survival-related data objects.

parseCoxph <- function(surv.formula, data){
  survdata <- data[!duplicated(data[, 'id']), ]; n <- nrow(survdata)
  ph <- coxph(surv.formula, survdata, x = T)
  
  survdata <- data.frame(id = survdata$id, ph$x, survtime = ph$y[, 1], status = ph$y[, 2])
  pS <- ncol(ph$x)
  sf <- survfit(ph)
  ft <- sf$time[sf$n.event >= 1] # failure times
  nev <- sf$n.event[sf$n.event >= 1] # Number of failures per failure time
  
  S <- lapply(1:n, function(i) ph$x[i, ])
  Delta <- as.list(survdata$status)
  
  #' #' Obtain form of S through call to ph ----  (If coxph object is entered instead of a formula)
  #' call <- as.character(ph$call)[2]
  #' S.covariates <- el(strsplit(call, '\\~'))[2]
  #' 
  #' #' Obtain survival object S, along with other information we want to know ----
  #' n <- length(ph$y)
  #' S <- cbind(data$id, model.matrix(as.formula(paste0('~', S.covariates, '-1')), data))
  #' S <- S[!duplicated.matrix(S), ]
  #' survtime <- unname(ph$y[,1])
  #' status <- unname(ph$y[,2])
  
  #' Cast certain items to List ----
  # Delta <- as.list(status)
  # S <- lapply(1:n, function(i) matrix(S[i, -1], nr = 1))
  
  #' Return ----
  list(
   survdata = survdata, ph = ph, n = n, S = S, ft = ft, nev = nev, Delta = Delta
  )
}

# surv.mod takes the fitted object and survival data from above function along with l0.init
# (from inits.R) and the (random effect) formulas, which help inform construction of F_x objects.
surv.mod <- function(ph, survdata, formulas, l0.init){
  uids <- unique(survdata$id); n <- length(uids)
  #' Unpack ph, data, formulas ----
  random <- as.character(formulas$random)[2]
  q <- ifelse(random == '1', 1, length(random) + 1)
  l0 <- l0.init
  
  #' Define two functions to construct F_i and F_u ----
  .getFi <- function(time, q = q){
    Fi <- matrix(NA, nr = 1, nc = q)
    for(i in 1:q) Fi[, i] <- time^(i - 1)
    Fi
  }
  .getFu <- function(times, q = q){
    sapply(1:q, function(i) times ^ (i - 1))
  }
  
  # Failure times
  ft <- coxph.detail(ph)$time
  
  # initialise empty stores
  l0i <- vector('numeric', n)
  Fu <- l0u <- surv.times <- vector('list', n)
  
  # get F_i
  Fi <- lapply(survdata$survtime, .getFi, q)
  
  #' loop over subjects -----
  for(i in uids){
    # slice this id
    survtime <- survdata[survdata$id == i, 'survtime']
    status <- survdata[survdata$id == i, 'status']
    # INDICES of survived times
    surv.times[[i]] <- which(ft <= survtime)
    # Fu, and l0u
    st <- ft[which(ft <= survtime)] # survived times
    if(length(st) > 0){
      Fu[[i]] <- .getFu(st, q)
      l0u[[i]] <- l0[which(ft <= survtime)]
    }else{ # Those who are censored before first failure time
      l0u[[i]] <- 0
      Fu[[i]] <- do.call(cbind, replicate(q, 0, simplify = F))
    }
    if(status == 1) l0i[i] <- l0[which(ft == survtime)] else l0i[i] <- 0
  }
  
  #' Return list ----
  return(list(
    ft = ft, nev = coxph.detail(ph)$nevent, surv.times = surv.times,
    l0 = l0, l0i = as.list(l0i), l0u = l0u, Fi = Fi, Fu = Fu
  ))
}


# Family-related ----------------------------------------------------------

setfamily <- function(family){
  if("function"%in%class(family)) family <- family()$family # step to ensure non-quoted arguments don't throw error.
  family <- match.arg(family, c('gaussian', 'binomial', 'poisson', 'negative.binomial'), several.ok = F)
  
  #' Switch statement
  switch(family,  
  gaussian = {
    log.likelihood <- gaussian_ll
    score.eta <- Score_eta_gauss
    sigma.update <- vare_update
    mu <- identity
  },
  binomial = {
    log.likelihood <- binomial_ll
    score.eta <- Score_eta_binom
    sigma.update <- NULL
    mu <- plogis
  },
  poisson = {
    log.likelihood <- poisson_ll
    score.eta <- Score_eta_poiss
    sigma.update <- NULL
    mu <- exp
  },
  negative.binomial = {
    log.likelihood <- negbin_ll
    score.eta <- Score_eta_negbin
    sigma.update <- theta_update
    mu <- exp
  })
  
  # Return a list of family-specific arguments
  return(list(log.likelihood = log.likelihood,
              score.eta = score.eta,
              sigma.update <- sigma.update,
              mu = mu
         ))
}


