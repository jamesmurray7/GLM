#' ###
#' Functions to compare my PBC (multivariate) fits with other software.
#' ###


.to3dp <- function(x) round(x, 3)
my.summary <- function(myfit){
  if(is.null(myfit$SE)) stop('Need to run EM with post.process = T')
  qz <- qnorm(.975)
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
  cat(paste0('Random effects variance-covariance matrix: \n'))
  print(.to3dp(D))
  cat('\n')
  
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

# Non-Gaussians -----------------------------------------------------------
library(JMbayes2)
JMbayes2.fit <- function(data, family, long.formula, random.formula, surv.formula, ...){
  survdata <- joineR::UniqueVariables(data, c('drug', 'survtime', 'status'), 'id')
  if("function"%in%class(family) & family()$family == 'gaussian'){
    M <- lme(long.formula, random = random.formula, data = data)
  }else{
    M <- mixed_model(long.formula, random.formula, data = data, family = family)
  }
  S <- coxph(surv.formula, data = survdata)
  jm(S, M, time_var = 'time', control = list(
    cores = 1, n_chains = 2
  ))
}


# Plotting functionality --------------------------------------------------
library(ggplot2)
plot.long <- function(data, Y, by){
  ggplot(data = data, aes(x = time, y = {{Y}}, group = id)) + 
    geom_line(alpha = .20) + 
    geom_smooth(aes(group = {{by}}, colour = {{by}}),
                method = 'lm', formula = y~splines::ns(x, df = 3)) +
    theme_light()
}

