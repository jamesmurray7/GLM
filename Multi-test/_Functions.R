#' ######
#' _Functions.R
#'   Collection of functions to help parse formula inputs and set-up data objects.
#' ######

# Parsing input formula ---------------------------------------------------
# (for longitudinal part), initial conditions found in inits.R
.parseFormula <- function(formula){
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

parseFormula <- function(formula){ # Wittled-down version which uses glmmTMB:::splitForm (a much more robust version of the above!)
  split <- glmmTMB:::splitForm(formula, allowFixedOnly = F)
  fixed <- split$fixedFormula
  random <- el(split$reTrmFormulas)
  
  #' Parse fixed effects
  response <- as.character(fixed)[2]
  fixed <- as.character(fixed)[3]
  if(grepl('splines\\:\\:|ns\\(|bs\\(', fixed)){
    attr(fixed, 'special') <- 'spline'
    if(grepl('ns\\(', fixed)) attr(fixed, 'spline type') <- 'natural'
    else if(grepl('bs\\(', fixed)) attr(fixed, 'spline type') <- 'basis'
    else stop('Unknown spline type')
  }else{
    attr(fixed, 'special') <- 'none'
  }
  # Flag splines as will likely be useful late wrhen setting up F_u, F_i, ftmat etc.^^^
  
  #' Parse random effects
  random.by <- as.character(random)[3]
  random <- as.character(random)[2]
  if(grepl('splines\\:\\:|ns\\(|bs\\(', random)){
    attr(random, 'special') <- 'spline'
    if(grepl('ns\\(', random)) attr(random, 'spline type') <- 'natural'
    else if(grepl('bs\\(', random)) attr(random, 'spline type') <- 'basis'
    else stop('Unknown spline type')
  }else{
    attr(random, 'special') <- 'none'
  }
  # Flag splines as will likely be useful later when setting up F_u, F_i, ftmat etc. ^^^
  
  return(list(
    response = response,
    fixed = fixed,
    random = random,
    random.by = random.by
  ))
}

.DataWorkhorse <- function(data, subj.id, fixed, random, response, what = 'X'){
  fixed.formula <- as.formula(paste0('~', fixed))
  random.formula <- as.formula(paste0('~', random))
  if(what == 'Y') rtn <- matrix(data[data$id == subj.id, response],nc = 1)
  else{
    if(!is.null(attr(fixed, 'special')) & (attr(fixed, 'special') == 'spline' | attr(random, 'special') == 'spline')){
      newData <- as.data.frame(cbind(id = data$id, model.matrix(fixed.formula, data)))
      newDataRandom <- as.data.frame(cbind(id = data$id, model.matrix(random.formula, data)))
      if(what == 'X') rtn <- as.matrix(newData[newData$id == subj.id, -1, drop = F])
      else if(what == 'Z') rtn <- as.matrix(newDataRandom[newDataRandom$id == subj.id, -1, drop = F])
    }else{ #' Otherwise proceed as usual.
      i.dat <- subset(data, id == subj.id)
      if(what == 'X') rtn <- model.matrix(fixed.formula, i.dat)
      else if(what == 'Z') rtn <- model.matrix(random.formula, i.dat)
    }
  }
  rtn
}

# Longitudinal objects ----------------------------------------------------
createDataMatrices <- function(data, formulas){
  K <- length(formulas)
  n <- length(unique(data$id))
  X <- Y <- Z <- vector('list', n)
  for(i in 1:n){
   ii <<- i
   X[[i]] <- lapply(formulas, function(f){
     .DataWorkhorse(data, subj.id = ii, fixed = f$fixed, random = f$random, response = f$response, what = 'X')
   }) 
   Z[[i]] <- lapply(formulas, function(f){
     .DataWorkhorse(data, subj.id = ii, fixed = f$fixed, random = f$random, response = f$response, what = 'Z')
   }) 
   Y[[i]] <- lapply(formulas, function(f){
     .DataWorkhorse(data, subj.id = ii, fixed = f$fixed, random = f$random, response = f$response, what = 'Y')
   }) 
  }
  
  list(X = X, Y = Y, Z = Z)
}

# Survival objects --------------------------------------------------------
# Parsing the survival formula and constructing all survival-related data objects.

parseCoxph <- function(surv.formula, data){
  survdata <- data[!duplicated(data[, 'id']), ]; n <- nrow(survdata)
  ph <- coxph(surv.formula, survdata, x = T)
  
  survdata <- data.frame(id = survdata$id, ph$x, survtime = ph$y[, 1], status = ph$y[, 2])
  pS <- ncol(ph$x)
  sf <- survfit(ph)
  ft <- sf$time[sf$n.event >= 1] # failure times
  nev <- sf$n.event[sf$n.event >= 1] # Number of failures per failure time
  
  Delta <- as.list(survdata$status)
  
  #' Return ----
  list(
   survdata = survdata, ph = ph, n = n, Delta = Delta
  )
}

# surv.mod takes the fitted object and survival data from above function along with l0.init
# (from inits.R) and the (random effect) formulas, which help inform construction of F_x objects.
surv.mod <- function(ph, survdata, formulas, l0.init){
  uids <- unique(survdata$id); n <- length(uids); K <- length(formulas)
  l0 <- l0.init; splines <- F
  
  # initialise empty stores
  l0i <- vector('numeric', n)
  l0u <- surv.times <- vector('list', n)
  
  # Failure times
  ft <- coxph.detail(ph)$time
  
  # First loop over subjects to create K-invariant objects.
  for(i in as.numeric(uids)){
    # slice this id
    survtime <- survdata[survdata$id == i, 'survtime']
    status <- survdata[survdata$id == i, 'status']
    # INDICES of survived times
    surv.times[[i]] <- which(ft <= survtime)
    # Fu, and l0u
    st <- ft[which(ft <= survtime)] # survived times
    if(length(st) > 0){
      l0u[[i]] <- l0[which(ft <= survtime)]
    }else{ # Those who are censored before first failure time
      l0u[[i]] <- 0
    }
    if(status == 1) l0i[i] <- l0[which(ft == survtime)] else l0i[i] <- 0
  }
  
  #' Second loop over formulas and subjects to create i x K object. ----
  
  sv <- lapply(formulas, function(formulas){
    if((!is.null(attr(formulas$random, 'special'))) & attr(formulas$random, 'special') == 'spline'){
      splines <- T
      survdata$time <- unname(ph$y[,1])
      newSurvdata <- as.data.frame(cbind(id = survdata$id, model.matrix(as.formula(paste0('~', formulas$random)), survdata)))
      q <- ncol(newSurvdata) - 1
      .getFi <- function(time, q = q){
        id <- which(survdata$time == time)[1] # Take the first one in case of ties -- same design matrix regardless
        as.matrix(newSurvdata[id, -1, drop = F])
      } 
      .getFu <- function(times, q = q){
        as.matrix(newSurvdata[match(times, survdata$time), -1, drop = F])
      }
      Fi <- lapply(survdata$survtime, .getFi)
      ftmat <- .getFu(ft)
    }else{ #' Case 2:: Anything else // Maybe more specials in future(?)
      random <- formulas$random
      q <- length(el(strsplit(random, '\\+|\\*|\\:')))
      #' Define two functions to construct F_i and F_u ----
      .getFi <- function(time, q = q){
        Fi <- matrix(NA, nr = 1, nc = q)
        for(i in 1:q) Fi[, i] <- time^(i - 1)
        Fi
      }
      .getFu <- function(times, q = q){
        out <- sapply(1:q, function(i) times ^ (i - 1))
        if(!"matrix"%in%class(out)) out <- t(out)
        out
      }
      # Generate F_i and ftmat
      Fi <- lapply(survdata$survtime, .getFi, q)
      ftmat <- .getFu(ft, q)
    }
    
    #' loop over subjects -----
    Fu <- vector('list', n)
    for(i in as.numeric(uids)){
      # slice this id
      survtime <- survdata[survdata$id == i, 'survtime']
      status <- survdata[survdata$id == i, 'status']
      # INDICES of survived times
      surv.times[[i]] <- which(ft <= survtime)
      # Fu, and l0u
      st <- ft[which(ft <= survtime)] # survived times
      if(length(st) > 0){
        Fu[[i]] <- .getFu(st, q)
      }else{ # Those who are censored before first failure time
        Fu[[i]] <- do.call(cbind, replicate(q, 0, simplify = F))
      }
    }
    list(Fu = Fu, Fi = Fi, ftmat = ftmat)
  })
  
  # Collate Fi, Fu, ftmat.
  Fi <- lapply(lapply(1:n, function(i){
    lapply(1:K, function(k){
      sv[[k]]$Fi[[i]]
    })
  }), function(x) do.call(cbind, x))
  
  Fu <- lapply(lapply(1:n, function(i){
    lapply(1:K, function(k){
      sv[[k]]$Fu[[i]]
    })
  }), function(x) do.call(cbind, x))
  
  ftmat <- lapply(1:K, function(k){
    sv[[k]]$ftmat
  })
  
  # Design matrix SS and rowvec S
  S <- lapply(1:n, function(i) ph$x[i, , drop = F])
  SS <- lapply(1:n, function(i){
    out <- apply(S[[i]], 2, rep, nrow(Fu[[i]]))
    if("numeric"%in%class(out)) out <- t(out)
    out
  })
  
  #' Return list ----
  return(list(
    ft = ft, ft.mat = ftmat, nev = coxph.detail(ph)$nevent, surv.times = surv.times,
    l0 = l0, l0i = as.list(l0i), l0u = l0u, Fi = Fi, Fu = Fu, Tis = survdata$survtime,
    S = S, SS = SS, q = q
  ))
}

difference <- function(params.old, params.new, type){
  if(type == 'absolute'){
    rtn <- max(abs(params.new - params.old))
  }else if(type == 'relative'){
    rtn <- max(
      abs(params.new - params.old)/(abs(params.old) + 1e-3)
    )
  }else{
    rtn <- NA
  }
  rtn
}

bind.bs<- function(bsplit){
  qmax <- max(sapply(bsplit, length)) # Maximum number of REs across all longitudinal responses.
  # Pad with zeros until each row is of length qmax
  step <- lapply(bsplit, function(b){
    l <- length(b)
    if(l<qmax) b <- c(b, rep(0, (qmax-l)))
    b
  })
  step <- do.call(rbind, step); colnames(step) <- NULL
  as.matrix(step)
}

my.summary <- function(myfit, printD = F){
  if(is.null(myfit$SE)) stop('Need to run EM with post.process = T')
  qz <- qnorm(.975)
  .to3dp <- function(x) round(x, 3)
  # Model fit info
  K <- length(myfit$ResponseInfo)
  responses <- lapply(sapply(myfit$ResponseInfo, strsplit, '\\s\\('), el, 1)
  families <- unlist(myfit$family)
  # Standard errors and parameter estimates.
  SE <- myfit$SE
  D <- myfit$co$D
  betas <- myfit$co$beta
  sigmas <- unlist(myfit$co$sigma)
  gammas <- myfit$co$gamma
  zetas <- myfit$co$zeta
  
  #' Random effects matrix
  if(printD){
    cat(paste0('Random effects variance-covariance matrix: \n'))
    print(.to3dp(D))
    cat('\n')
  }
  
  # Longitudinal parts
  MakeTables <- lapply(1:K, function(k){
    
    beta <- betas[grepl(responses[[k]], names(betas))]
    sigma <- setNames(sigmas[k], paste0(responses[[k]], '_var.e'))
    
    cat(paste0(responses[[k]], ' (', families[k], '): \n'))
    my <- c(beta)
    if(sigma != 0.0) my <- c(my, sigma)
    parameter <- names(my)
    
    rSE <- SE[match(names(my), names(SE))]#SE associated with these coeffs
    
    lb <- my - qz * rSE; ub <- my + qz * rSE
    
    z <- my/rSE
    p <- 2 * (pnorm(abs(z), lower.tail = F))
    
    this.out <- setNames(data.frame(.to3dp(my), .to3dp(rSE), .to3dp(lb), .to3dp(ub), round(p, 3)),
                         c('Estimate', 'SE', '2.5%', '97.5%', 'p-value'))
    print(this.out)
    cat('\n')
  })
  
  #' Survival
  cat('Event-time sub-model: \n')
  # Rename gammas?
  survs <- c(zetas, gammas)
  surv.SE <- SE[match(names(survs), names(SE))]
  new.gammas <- sapply(1:K, function(k) gsub('\\_\\d?\\d', paste0('_', unlist(responses)[k]), names(gammas)[k]))
  names(survs)[grepl('gamma', names(survs))] <- new.gammas
  
  lb <- survs - qz * surv.SE; ub <- survs + qz * surv.SE
  z <- survs/surv.SE
  p <- 2 * (pnorm(abs(z), lower.tail = F))
  
  surv.out <- setNames(data.frame(.to3dp(survs), .to3dp(surv.SE), .to3dp(lb), .to3dp(ub), round(p, 3)),
                       c('Estimate', 'SE', '2.5%', '97.5%', 'p-value'))
  print(surv.out)
  
  cat(paste0('\nApproximate EM algorithm took ', round(myfit$EMtime, 2), ' seconds and total computation time was ', round(myfit$totaltime, 2), ' seconds.\n'))
  invisible(1+1)
}

fit2xtab <- function(fit, max.row = NULL){
  if(is.null(fit$SE)) stop('Need to run EM with post.process = T')
  qz <- qnorm(.975)
  .to3dp <- function(x) round(x, 3)
  
  # Model fit info
  K <- length(fit$ResponseInfo)
  responses <- lapply(sapply(fit$ResponseInfo, strsplit, '\\s\\('), el, 1)
  families <- unlist(fit$family)
  # Standard errors and parameter estimates.
  SE <- fit$SE
  D <- fit$co$D
  betas <- fit$co$beta
  sigmas <- unlist(fit$co$sigma)
  gammas <- fit$co$gamma
  zetas <- setNames(fit$co$zeta, paste0('zeta_', 1:length(fit$co$zeta)))
  
  MakeTables <- lapply(1:K, function(k){
    
    nb <- names(betas)[grepl(responses[[k]], names(betas))]
    nb2 <- paste0('beta_{', k, (seq(0, (length(nb) - 1))), '}')
    beta <- setNames(betas[grepl(responses[[k]], names(betas))], nb2)
    if(sigmas[k] != 0) sigma <- setNames(sigmas[k], paste0('sigma^2_', k)) else sigma <- NULL
    gamma <- setNames(gammas[k], paste0('gamma_', k))
    
    kk <- c(beta, sigma, gamma)
    if(sigmas[k] != 0) sigma.name.lookup <- paste0(responses[[k]], '_var.e') else sigma.name.lookup <- NULL
    kk.names.lookup <- c(nb, sigma.name.lookup, names(gamma))
    
    kSE <- SE[match(kk.names.lookup, names(SE))]#SE associated with these coeffs
    
    lb <- kk - qz * kSE; ub <- kk + qz * kSE
    
    this.out <- setNames(data.frame(.to3dp(kk), .to3dp(kSE), .to3dp(lb), .to3dp(ub)),
                         c('Estimate', 'SE', '2.5%', '97.5%'))
    
    this.out
  })
  
  tab <- do.call(rbind, MakeTables)
  # Append zeta terms to bottom, time invariant so report separately
  SEz <- SE[grepl('^zeta', names(SE))]; lb <- zetas - qz * SEz; ub <- zetas + qz * SEz
  zet <- data.frame(.to3dp(zetas), .to3dp(SEz), .to3dp(lb), .to3dp(ub))
  names(zet) <- names(tab)
  tab <- rbind(tab, zet)
  nr <- nrow(tab)
  
  tab$Parameter <- rownames(tab)
  tab$Parameter <- paste0('$\\', tab$Parameter, '$')
  
  tab2 <- as.data.frame(cbind(Parameter = tab$Parameter, apply(tab[, -5], 2, function(x) format(round(x, 3), nsmall = 3))))
  tab3 <- cbind(Parameter = tab$Parameter, `Mean (SE)` = paste0(tab2$Estimate, ' (', tab2$SE, ')'), 
                `95% CI` = paste0('[', tab2$`2.5%`, ', ', tab2$`97.5%`, ']'))
  
  # Splitting out into multiple columns -->
  if(nr > 15 && is.null(max.row)){
    cat('Consider breaking at a certain number of rows and presenting a "wider" table.\nnrows: ', nr, '\n') 
  }
  
  if(!is.null(max.row)){
    if(nr <= max.row) stop('max.row must exceed the number of rows in output table: ', nr, '.\n')
    
    # Work out how many 'cbinds' we'll need to do.
    num.splits <- nr %/% max.row
    nr.each <- suppressWarnings(split(1:nr, rep(1:num.splits, each = ceiling(nr/num.splits))))
    split.tab <- lapply(nr.each, function(k){
      x <- tab3[k,]
      while(nrow(x) < max.row){
        x <- rbind(x, c('-','-','-'))
      }
      x
    })
    
    tab3 <- do.call(cbind, split.tab)
  }
  
  xt <- xtable::xtable(tab3,
                       caption = paste0("Elapsed time for approximate EM algorithm to converge and SE calculation was ", round(fit$EMtime + fit$postprocess.time, 2), " seconds."))
  
  print(xt,
        include.rownames = FALSE,
        sanitize.text.function = identity)
  
}

log.lik <- function(coeffs, dmats, b, surv, sv, l0u, l0i, gamma.rep, beta.inds, b.inds, K, q, family){
  # joint density - REs; Marginal ll of Y(s) added to f(T_i,\Delta_i|...).
  Y <- dmats$Y; X <- dmats$X; Z <- dmats$Z
  beta <- coeffs$beta; D <- coeffs$D; sigma <- coeffs$sigma; zeta <- coeffs$zeta
  S <- sv$S; SS <- sv$SS; Delta <- surv$Delta
  Fu <- sv$Fu; Fi <- sv$Fi; 

  #' "Full" joint density ----
  ll1 <- mapply(function(Y, X, Z, b, S, SS, Fu, Fi, Delta, l0i, l0u){
    joint_density(b = b, Y = Y, X = X, Z = Z, beta = beta, D = D, sigma = sigma, family = family, 
                  Delta = Delta, S = S, Fi = Fi, l0i = l0i, SS = SS, Fu = Fu, haz = l0u, gamma_rep = gamma.rep, zeta = zeta,
                  beta_inds = beta.inds, b_inds = b.inds, K = K) * -1
  }, Y = Y, X = X, Z = Z, b = b, S = S, SS = SS, Fu = Fu, Fi = Fi, Delta = Delta, l0i = l0i, l0u = l0u)
  
  #' Only the joint density associated with REs ----
  ll2 <- mapply(function(b) dmvnrm_arma_fast(t(b), rep(0, q), D), b = b)
  #' Calculate the marginal ll of the Yks and the survival density.
  out <- sum(ll1 - ll2)
  
  #' Calculate AIC and BIC
  N <- sum(colSums(do.call(rbind, lapply(Y, function(y) sapply(y, length)))))
  P <- sum(sapply(1:K, function(k) ncol(X[[1]][[k]]) - 1))
  Ps <- ncol(S[[1]])
  df <- P + Ps + 2 * K + (q * (q + 1)/2)
  
  attr(out, 'df') <- df; attr(out, 'N') <- N
  attr(out, 'AIC') <- -2 * c(out) + 2 * df
  attr(out, 'BIC') <- -2 * c(out) + log(N) * df
  out
}
