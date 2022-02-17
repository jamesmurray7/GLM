#' ###
#' Printing a summary for a given EM fit
#' ###

summary.mix <- function(fit, ci = 95){
  if(is.null(fit$SE)) stop('Need to run EM with post.process = T')
  SE <- fit$SE
  NB.ind <- any(grepl('^NB\\_', names(SE)))
  coeff <- fit$coeffs
  names.to.report <- !grepl('^D', names(SE)) # report everything but covariance matrix
  
  # Establishing quantile stuff
  alpha <- (100 - ci)/100
  qz <- qnorm(1-alpha/2)
  
  lower.name <- paste0(100 * alpha/2, '%')
  upper.name <- paste0(100 * (1-alpha/2), '%')
  
  # Longitudinal part
  coeff.long <- c(c(coeff$beta), coeff$var.e)
  if(NB.ind) coeff.long <- c(coeff.long, coeff$theta)
  SE <- SE[names.to.report]
  longit.ind <- grepl('^G\\_|^B\\_|^P\\_|^NB\\_|var.e', names(SE))
  SE.long <- SE[longit.ind]
  names(coeff.long) <- names(SE.long)
  
  lb <- coeff.long - qz * SE.long
  ub <- coeff.long + qz * SE.long
  
  z <- coeff.long/SE.long
  p <- 2 * (pnorm(abs(z), lower.tail = F))
  
  longpart <- setNames(
    data.frame(coeff.long, SE.long, lb, ub, z, p),
    c('Coefficient', 'SE', lower.name, upper.name, 'z-value', 'p-value')
  )
  
  cat(print(longpart))
  
  # Survival Part
  
  
}
