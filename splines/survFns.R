#' #######
#' Functions to obtain data objects from the survival sub-model,
#' given a coxph fit, data and initial conditions for the baseline hazard (optional).
#' Allows the user to specify degree of basis used here.
#' ---
#' These are:
#' GLOBAL:::
#' l0: hazard associated with failure times
#' ft: Failure times (non-censored)
#' SUBJECT-SPECIFIC:::
#' Fu: design matrix on all failure times survived by subject i
#' Fi: design matrix of subject i's failure time (1, Ti).
#' l0u: hazard UP TO subject i's failure time
#' l0i: hazard associated with subject i's failure time
#' ---
#' #######

surv.mod <- function(cph, data, l0.init = NULL, degree = 3){
  uids <- unique(data$id)
  if("coxph.null"%in%class(cph)) message("Null model")
  # Survfit
  sf <- summary(survfit(cph))
  
  # initialise empty stores
  Fi <- matrix(NA, nr = length(uids), nc = (degree + 1))
  Di <- l0i <- numeric(length(uids))
  Fu <- l0u <- surv.times <- vector('list', length(uids))
  
  if(is.null(l0.init)){
    l0 <- diff(c(0, sf$cumhaz))
  }else{
    l0 <- l0.init
  }
  
  # Define basis
  .Tis <- unique(data[, 'survtime'])
  .ft <- unique(data[data$status == 1, 'survtime'])
  survbasis <- bs(.Tis, degree = degree)
  
  # loop
  for(i in uids){
    i.dat <- unique(subset(data, id == i, c('survtime', 'status')))
    survtime <- i.dat$survtime; status <- i.dat$status
    if(nrow(i.dat) != 1) stop('Error, non-unique failure times reported for subject ', i, '.')
    # Failure indicator
    Di[i] <- status
    # Vector of time (indices) survived
    surv.times[[i]] <- which(sf$time <= survtime)
    # This subject's Fi: (1, Ti) structure
    Fi[i,] <- c(1, survbasis[which(.Tis == survtime), ])
    # This subject's design matrix Fu on all (ordered failure) times survived.
    Fu[[i]] <- cbind(1, as.matrix(survbasis[match(.ft[which(.ft <= survtime)], .Tis), ]))
    # This subject's hazard vector
    l0u[[i]] <- l0[which(sf$time <= survtime)]
    # The individual hazard \lambda(Ti) at failure time; zero if censored
    if(status == 1) l0i[i] <- l0[which(sf$time == survtime)] else l0i[i] <- 0
    # Check if censored before first failure time. Set \lambda(u) to zero and set Fu to be empty matrix of appropriate dimension
    if(status == 0 & survtime <= min(sf$time)){ l0u[[i]] <- 0; Fu[[i]] <- cbind(1, matrix(0, nr = 1, nc = degree)) }
  }
  
  nev <- c(); surv.ids <- list()
  p <- 1
  for(i in sf$time){
    nev[p] <- length(unique(data[which(data$survtime == i),]$id))
    surv.ids[[p]] <- unique(data[which(data$survtime >= i),]$id)
    p <- p+1
  }
  
  # output
  return(list(
    ft = sf$time,
    l0 = l0,
    nev = nev,
    surv.ids = surv.ids,
    surv.times = surv.times,
    l0i = l0i,
    Di = Di,
    Fi = Fi,
    Fu = Fu,
    l0u = l0u,
    basis = structure(cbind(.Tis, survbasis),
                      dimnames = list(as.character(1:nrow(survbasis)),
                                      c('Ti', paste0('basis', 1:degree)))
                      )
  ))
}
