library(tidyverse)
theme_csda <- function(base_family = "Arial", base_size = 12){
  theme_light(base_family = base_family, base_size = base_size) %+replace%
    theme(
      strip.background = element_blank(),
      strip.text = element_text(size = 11, colour = "black", vjust = 1),
      axis.text = element_text(size = 9),
      axis.title = element_text(size = 11),
      legend.background = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
}
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
  f <- fit$coefs.long[, 1:2]
  s <- rbind(fit$coefs.surv[-1, ], fit$coefs.surv[1, ])[, 1:2] # swap order gamma zeta
  row.names(s) <- rev(row.names(fit$coefs.surv))
  var <- cbind(Value = fit$coefs.sigma, `Std.Err` = 0)
  
  data.frame(rbind(f, s, var^2), method = 'joineRML')
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
    jML <- jML %>% select(method, param, estimate = Value)
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


# Tabulation --------------------------------------------------------------
toXdp <- function(x, X) format(round(x, X), nsmall = X)
tabulate <- function(fit, type, f, K){ # fit a sub-list
  type <- match.arg(type, c('aEM', 'jML', 'JMb', 'INLA'), several.ok = F)
  f <- match.arg(f, c('gaussian', 'poisson', 'binomial'), several.ok = F)
  
  targets <- get.targets(f, K)
  
  if(!is.null(fit)){ # Parse
    if(type == 'aEM') fit <- as.data.frame(parse.aEM(fit))
    if(type == 'jML') fit <- as.data.frame(parse.jML(fit))
    if(type == 'JMb') fit <- as.data.frame(parse.JMb(fit))
    if(type == 'INLA') fit <- as.data.frame(parse.INLA(fit))
  }else return(NULL)
  fit$target <- targets
  # Specific parsing...
  r <- row.names(fit); qz <- qnorm(.975)
  # cat(r,'\n',names(targets),'\n\n')
  switch(type,
         'aEM' = {
           fit$is.cp <- fit$lb <= targets & fit$ub >= targets
           fit$bias <- fit$Omega - fit$target
           fit$param <- names(targets)
           fit <- fit %>% select(param = param, target = target, estimate = Omega, SE, bias, is.cp)
         },
         'jML' = {
           lb <- fit$Value - qz * fit$`Std.Err`
           ub <- fit$Value + qz * fit$`Std.Err`
           fit$is.cp <- lb <= targets & ub >= targets
           fit$bias <- fit$Value - fit$target
           fit$param <- names(targets)
           fit <- fit %>% select(param = param, target = target, estimate = Value, SE = `Std.Err`, bias, is.cp)
         },
         'JMb' = {
           fit$is.cp <- fit$`2.5%` <= targets & fit$`97.5%` >= targets
           fit$bias <- fit$Mean - fit$target
           fit$param <- names(targets)
           fit <- fit %>% select(param = param, target = target, estimate = Mean, SE = StDev, bias, is.cp)
         },
         'INLA' = {
           fit$is.cp <- fit$`0.025quant` <= targets & fit$`0.975quant` >= targets
           fit$bias <- fit$mean - fit$target
           fit$param <- names(targets)
           fit <- fit %>% select(param = param, target = target, estimate = mean, SE = sd, bias, is.cp)
         })
  return(fit)
}
 
tabulate_wrapper <- function(fit, type, f, K){ # fit a list of lists
  tab <- do.call(rbind, lapply(fit, tabulate, type = type, f = f, K = K))
  # print(head(tab))
  if(type!='jML'){
    out1 <- tab %>% 
      group_by(param) %>% 
      summarise(target = unique(target),
                mean = mean(estimate),
                sd = mean(SE),
                bias.m = mean(bias),
                bias.s = sd(bias),
                CP = sum(is.cp)/n(),
                conv = n()/100) %>% 
      ungroup
  }else{
    out1 <- tab %>% 
      filter(!grepl('^sigma', param)) %>% 
      group_by(param) %>% 
      summarise(target = unique(target),
                mean = mean(estimate),
                sd = mean(SE),
                bias.m = mean(bias),
                bias.s = sd(bias),
                CP = sum(is.cp)/n(),
                conv = n()) %>% 
      ungroup
    # print(out1)
    sigmas <- tab %>% 
      filter(grepl('^sigma', param)) %>% 
      group_by(param) %>% 
      mutate(s = sd(estimate)) %>% 
      ungroup() %>% 
      mutate(is.cp = (estimate - qnorm(.975) * s) <= target & 
               (estimate + qnorm(.975) * s) >= target) %>% 
      group_by(param) %>% 
      summarise(
        target = unique(target),
        mean = mean(estimate),
        sd = s,
        bias.m = mean(bias),
        bias.s = sd(bias),
        CP = sum(is.cp)/n(),
        conv = n()/100
      ) %>% 
      ungroup %>% distinct
    # print(sigmas)
    out1 <- rbind(out1, sigmas) %>% arrange(param)
  }
  
  out1 %>% mutate_at(c('target', 'mean', 'sd', 'bias.m', 'bias.s'), ~ toXdp(.x, 3)) %>% 
    mutate_at(c('CP', 'conv'), ~ toXdp(.x, 2)) %>% 
    mutate(mean.SD = paste0(mean, ' (', sd, ')'),
           bias.SD = paste0(bias.m, ' (', bias.s, ')')) %>% 
    select(param, target, mean.SD, bias.SD, CP, conv) %>% 
    rename_at(c('mean.SD', 'bias.SD', 'CP', 'conv'),
              ~ paste0(type, '-', .x))
}
library(xtable)
tab.to.xtab <- function(tab, mean.col = F){
  tab$param <- gsub('\\[', '_{', tab$param)
  tab$param <- gsub('\\]', '}', tab$param)
  cn <- colnames(tab)
  convs <- tab[1,grepl('conv', cn)]
  print(convs)
  
  # Xtable multicolumn stuff
  if(any(grepl('jML', cn))){
    methods <- c('Approximate EM', 'joineRML', 'JMbayes2', 'INLAjoint')
  }else{
    methods <- c('Approximate EM', 'JMbayes2', 'INLAjoint')
  }
  addtorow <- list()
  addtorow$pos <- list(-1)
  addtorow$command <- paste0(paste0('& \\multicolumn{3}{c}{', methods, '}', collapse=''), '\\\\')
  
  if(mean.col){
    tab <- tab %>% 
      select(-ends_with('conv')) %>% 
      rename_at(vars(`aEM-mean.SD`:`aEM-CP`), ~gsub('aEM\\-','a',.x)) %>% 
      rename_at(vars(`JMb-mean.SD`:`JMb-CP`), ~gsub('JMb\\-','b',.x)) %>% 
      rename_at(vars(`INLA-mean.SD`:`INLA-CP`), ~gsub('INLA\\-','c',.x))
    if(any(grepl('jML', cn))) tab <- tab %>% rename_at(vars(`jML-mean.SD`:`jML-CP`), ~gsub('jML\\-','d',.x)) 
    addtorow$command <- paste0(paste0('& \\multicolumn{3}{c}{', methods, '}', collapse=''), '\\\\')
  }else{
    tab <- tab %>% select(-ends_with('mean.SD'))
    tab <- tab %>% 
      select(-ends_with('conv')) %>% 
      rename_at(vars(`aEM-bias.SD`:`aEM-CP`), ~gsub('aEM\\-','a',.x)) %>% 
      rename_at(vars(`JMb-bias.SD`:`JMb-CP`), ~gsub('JMb\\-','b',.x)) %>% 
      rename_at(vars(`INLA-bias.SD`:`INLA-CP`), ~gsub('INLA\\-','c',.x))
    if(any(grepl('jML', cn))) tab <- tab %>% rename_at(vars(`jML-bias.SD`:`jML-CP`), ~gsub('jML\\-','d',.x)) 
    addtorow$command <- paste0(paste0('& \\multicolumn{2}{c}{', methods, '}', collapse=''), '\\\\')
  }
  
  xt <- tab %>% 
    mutate(
      param = paste0('$\\', param, '=', target,'$')
    ) %>% select(-target) %>% 
    rename_at(vars(contains('bias.SD')), ~gsub('bias.SD', 'Bias (SD)', .x)) %>% 
    xtable::xtable()
    
  print(xt,
        include.rownames = FALSE,
        sanitize.text.function = identity,
        add.to.row = addtorow)
}
