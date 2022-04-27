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
