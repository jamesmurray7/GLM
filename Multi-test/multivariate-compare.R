#' ###
#' Functions to compare my PBC (univariate) fits with joineR, joineRML and JMbayes2.
#' ###

library(joineR)
joineR.fit <- function(data, Y, long.formula, surv.formula, long.covs, baseline.covs){
  survdata <- UniqueVariables(data, var.col = c('survtime', 'status'), id.col = 'id')
  longdata <- data[, c('id', Y, long.covs)]
  covdata <- UniqueVariables(data, baseline.covs, id.col = 'id')
  joint.data <- jointdata(longitudinal = longdata,
                          baseline = covdata,
                          survival = survdata,
                          id.col = 'id', time.col = 'time')
  fit <- joint(data = joint.data,
               long.formula = long.formula,
               surv.formula = surv.formula,
               model = 'intslope', sepassoc = F, gpt = 3)
  return(fit)
}

library(joineRML)
joineRML.fit <- function(data, long.formula, random.formula, surv.formula, ...){
  mjoint(
    formLongFixed = list('1' = long.formula),
    formLongRandom = list('1' = random.formula),
    formSurv = surv.formula,
    data = data,
    timeVar = 'time',
    control = list(convCrit = 'rel', type = 'sobol', tol2 = 1e-2, ...)
  )
}

# Compare Gaussians -------------------------------------------------------
compare.gaussians <- function(myfit, jRfit, jRMLfit, Y){
  #' Mine
  my <- c(vech(myfit$co$D), c(myfit$co$beta), myfit$co$gamma, myfit$co$zeta)
  if(myfit$co$sigma != 0) my <- c(my, myfit$co$sigma)
  parameter <- names(myfit$SE)
  SE <- myfit$SE
  lb <- my - qnorm(.975) * SE; ub <- my + qnorm(.975) * SE
  
  #' joineR point-estimates
  j.ests <- c(vech(jRfit$sigma.u), c(jRfit$coefficients$fixed$longitudinal)$b1, jRfit$coefficients$latent,
              jRfit$coefficients$fixed$survival, jRfit$sigma.z)
  
  my.out <- setNames(data.frame(parameter, my, SE, lb, ub, j.ests),
                     c('Parameter', 'Estimate', 'SE', '2.5%', '97.5%', 'joineR'))
  row.names(my.out) <- NULL
  #' joineR
  
  
  #' joineRML
  jML <- summary(jRMLfit)
  jML.out <- as.data.frame(rbind(jML$coefs.long, jML$coefs.surv))
  jML.out$`2.5%` <- jML.out$Value - qnorm(.975) * jML.out$Std.Err
  jML.out$`97.5%` <- jML.out$Value + qnorm(.975) * jML.out$Std.Err
  jML.out$Parameter <- rownames(jML.out)
  jML.out <- with(jML.out, setNames(data.frame(Parameter, Value, `Std.Err`, `2.5%`, `97.5%`),
                                    c('Parameter', 'Estimate', 'SE', '2.5%', '97.5%')))
  
  #' Print out
  cat(paste0('Model comparisons for Gaussian fits of ', Y, ':\n'))
  cat("Approximate EM Algorithm and joineR point-estimates:\n")
  print(my.out)
  cat('\nComputation times: \n')
  cat(paste0('EM took ', round(myfit$EMtime, 2), ' seconds and total computation time was ', round(myfit$totaltime, 2), ' seconds.\n'))
  cat('\njoineRML fit:\n')
  print(jML.out)
  cat('\nComputation times: \n')
  units <- attr(jML$comp.time, 'units')
  cat(paste0('EM took ', round(as.numeric(jML$comp.time[2]), 2), ' ', units, ' and total computation time was ', round(as.numeric(jML$comp.time[1]), 2), ' ', units,'.\n'))
  
  invisible(Y)
}

.to3dp <- function(x) round(x, 3)
my.summary <- function(myfit, Yname = NULL){
  if(!is.null(Yname)) cat('Approximate EM fit for ', Yname,'\n\n')
  if(is.null(myfit$SE)) stop('Need to run EM with post.process = T')
  SE <- myfit$SE
  
  # Establishing quantile stuff
  qz <- qnorm(.975)
  
  # Longitudinal part
  my <- c(vech(myfit$co$D), c(myfit$co$beta), myfit$co$gamma, myfit$co$zeta)
  if(myfit$co$sigma != 0) my <- c(my, myfit$co$sigma)
  parameter <- names(myfit$SE)
  lb <- my - qz * myfit$SE; ub <- my + qz * SE
  
  z <- my/SE
  p <- 2 * (pnorm(abs(z), lower.tail = F))
  
  my.out <- setNames(data.frame(parameter, .to3dp(my), .to3dp(SE), .to3dp(lb), .to3dp(ub), round(p, 3)),
                     c('Parameter', 'Estimate', 'SE', '2.5%', '97.5%', 'p-value'))
  row.names(my.out) <- NULL
  
  print(my.out)
  cat('\nComputation times: \n')
  cat(paste0('EM took ', round(myfit$EMtime, 2), ' seconds and total computation time was ', round(myfit$totaltime, 2), ' seconds.\n'))
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

