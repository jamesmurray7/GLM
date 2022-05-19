library(tidyverse)
# Functions for elapsed time ----------------------------------------------
elapsed.aEM <- function(fit, what = 'EM'){
  if(what == 'EM') return(fit$EMtime + fit$postprocess.time)
  else if(what == 'total') return(fit$totaltime)
}

elapsed.jML <- function(fit, what = 'EM'){
  if(what == 'EM'){
    x <- fit$comp.time[2]
    if(attr(x, 'units') == 'mins') x <- x * 60 else x <- x
    rtn <- as.numeric(x)
  }else if(what == 'total'){
    x <- fit$comp.time[1]
    if(attr(x, 'units') == 'mins') x <- x * 60
    rtn <- as.numeric(x)
  }
  rtn
}

elapsed.JMb <- function(fit){
  if(!is.null(fit)) return(fit$comp.time[3]) else return(NA)
}

elapsed.INLA <- function(fit){
  if(!is.null(fit)) return(unname(fit$comp.time[2])) else return(NA)
}

# Functions for object parsing --------------------------------------------

parse.aEM <- function(fit){
  co <- fit$coeffs
  beta <- co$beta
  if(any(unlist(co$sigma) != 0)) sig <- unlist(co$sigma) else sig <- NULL
  gam <- co$gamma
  ze <- co$zeta
  SE <- fit$SE[!grepl('^D', names(fit$SE))]
  Omega <- setNames(c(beta, gam, ze, sig), names(SE))
  
  lb <- Omega - qnorm(.975) * SE; ub <- Omega + qnorm(.975) * SE
  return(cbind(Omega, SE, lb, ub))
}

parse.JMb <- function(fit){
  o1 <- fit$Outcome1
  if(!is.null(fit$Outcome2)) o2 <- fit$Outcome2 else o2 <- NULL 
  if(!is.null(fit$Outcome3)) o3 <- fit$Outcome3 else o3 <- NULL
  o <- rbind(o1, o2, o3)
  r <- row.names(o)
  vars <- o[grepl('sigma', r),]
  o <- rbind(o[!grepl('sigma', r), ])
  s <- rbind(fit$survival[-1, ], fit$survival[1,]) # swap
  as.matrix(rbind(o, s, vars^2)[,c('Mean', 'StDev', '2.5%', '97.5%')])
}

parse.INLA <- function(fit){
  f <- fit$fixed
  # Ensure residual variance goes at the end
  rf <- row.names(f)
  vars <- f[grepl('\\(var\\)', rf),]
  f <- f[!grepl('\\(var\\)', rf),]
  z <- fit$survival[-c(1,2), ]
  g <- fit$gamma
  method <- 'INLA'
  cbind(rbind(f, g, z, vars), method)[,c('mean', 'sd', '0.025quant', '0.975quant', 'method')]
}

parse.jML <- function(fit){
  f <- fit$coefs.long[, 1]
  s <- rbind(fit$coefs.surv[-1, ], zeta = fit$coefs.surv[1, ])[, 1] # swap order gamma zeta
  var <- fit$coefs.sigma
  
  data.frame(est = c(f, s, var^2), method = 'joineRML')
}


# Obtaining targets for simulation study ----------------------------------
get.targets <- function(f, K){
  # Setting names
  if(K == 1){
    nm <- c(paste0('beta[', 0:3, ']'), 'gamma[1]', 'zeta') 
  }else{
    b <- paste0('beta[', 1:K)
    g <- paste0('gamma[', 1:K, ']')
    b <- paste0(rep(b, each = 4), 0:3, ']')
    nm <- c(b, g, 'zeta')
  }
  
  if(f == 'poisson'){
    targets <- setNames(c(do.call(c, replicate(K, c(2, -0.1, 0.1, -0.2), simplify = F)),
                          c(1, -1, 1)[1:K] * 0.3,
                          -0.2),
                        nm)
  }else if(f == 'binomial'){
    targets <- setNames(c(do.call(c, replicate(K, c(1, -0.1, 0.1, -0.2), simplify = F)),
                          c(1, -1, 1)[1:K] * 0.3,
                          -0.2),
                        nm)
  }else if(f == 'gaussian'){
    if(K > 1) nm <- c(nm, paste0('sigma[', 1:K, ']^2')) else nm <- c(nm, 'sigma^2')
    targets <- setNames(c(do.call(c, replicate(K, c(0.2, -0.1, 0.1, -0.2), simplify = F)),
                          c(1, -1, 1)[1:K] * 0.5,
                          -0.2, rep(0.16, K)),
                        nm)
  }else{
    stop("'f' must be gaussian/binomial/poisson")
  }
  
  targets
}


# Functions for plotting --------------------------------------------------
plot.ests <- function(my.in = NULL, jML.in = NULL, JMb.in = NULL, INLA.in = NULL,
                      f = 'gaussian', K){
  
  if(f != 'gaussian' & !is.null(jML.in)) stop("joineRML only Gaussian!")
  
  targets <- get.targets(f, K)
  print(targets)

  # Load in data if object is supplied
  if(!is.null(my.in)){
    my <- as.data.frame(do.call(rbind, lapply(my.in, parse.aEM)))
    my$method <- 'Approximate EM'
    my$param <- rep(names(targets), sum(!unlist(lapply(my.in, is.null))))
    my <- my %>% select(method, param, estimate = Omega)
  }else my <- NULL
  if(!is.null(jML.in)){
    jML <- as.data.frame(do.call(rbind, lapply(jML.in, function(x) if(!is.null(x)) parse.jML(x))))
    jML$method <- 'joineRML'
    jML$param <- rep(names(targets), sum(!unlist(lapply(jML.in, is.null))))
    jML <- jML %>% select(method, param, estimate = est)
  }else jML <- NULL
  if(!is.null(JMb.in)){
    JMb <- as.data.frame(do.call(rbind, lapply(JMb.in, function(x) if(!is.null(x)) parse.JMb(x))))
    JMb$method <- 'JMbayes2'
    JMb$param <- rep(names(targets), sum(!unlist(lapply(JMb.in, is.null))))
    JMb <- JMb %>% select(method, param, estimate = Mean)
  }else JMb <- NULL
  if(!is.null(INLA.in)){
    INLA <- as.data.frame(do.call(rbind, lapply(INLA.in, function(x) if(!is.null(x)) parse.INLA(x))))
    INLA$method <- 'INLA'
    INLA$param <- rep(names(targets), sum(!unlist(lapply(INLA.in, is.null))))
    INLA <- INLA %>% select(method, param, estimate = mean)
  }else INLA <- NULL

  ests <- rbind(my, jML, JMb, INLA) %>% 
    left_join(., data.frame(param = names(targets), target = targets), 'param')
  
  print(ests)
  
  ggplot(ests, aes(x = estimate, colour = method)) +
    geom_vline(aes(xintercept = target), lty = 5, colour = 'grey') +
    geom_density() + 
    facet_wrap(~ param, scales = 'free', labeller = label_parsed) + 
    labs(x = 'Estimate', y = NULL, colour = NULL,
         title = paste0('K = ', K, ' ', str_to_sentence(f))) +
    theme_light() + 
    theme(
      strip.background = element_blank(),
      strip.text = element_text(size = 12, colour = 'black')
    )
}


