#' #########
#' Getting handles dispersions of ADNI variables.
#' #########
rm(list=ls())
load('../../PBC-case-study/ADNI.RData') # Same across Linux/Mac(!)
source('EM.R')
library(dplyr)
# Define a function that 'resets' the data --------------------------------
newadni <- function(Y){
  y <- adni[,Y]
  inds.to.remove <- is.na(y) | y < 0
  if(length(inds.to.remove) > 0){
    rtn <- adni[!inds.to.remove, ]
    rtn <- rtn %>% 
      group_by(id) %>% 
      mutate(new_id = cur_group_id()) %>% 
      ungroup %>% 
      select(-id) %>% 
      rename(id = new_id) %>% 
      as.data.frame()
  }else{
    rtn <- adni
  }
  rtn
}

# Checking for over/under-dispersion in count responses -------------------
gen.summary <- function(response){
  cat(response, '----\n\n')
  d <- newadni(response)
  
  vars <- tapply(d[,response], d$id, var, na.rm = T)
  means <- tapply(d[,response], d$id, mean, na.rm = T)
  
  lengths <- median(with(d, tapply(time, id, length)))
  
  nas <- unname(which(is.na(vars) | vars == 0))
  vars <- vars[-nas]; means <- means[-nas]
  
  cat(sprintf('%d median profile length, mean/var summary:\n', lengths))
  print(summary(means/vars))
  
  cat('Making plot...\n\n')
  plot(vars~means,pch=20, main = response, ylab = 'Var[Y]', xlab = bquote(bar(X[i])))
  abline(0, 1, col = 'red', lwd = 0.75)
  
}

responses <- c('ADAS11', 'ADAS13', 'MMSE', 'RAVLT.immediate', 'RAVLT.learning', 'RAVLT.forgetting', 'FAQ')
sink('../ADNI/disps.txt')
invisible(sapply(responses, gen.summary))
sink()


# Function to find delta... -----------------------------------------------
disp.formula <- ~1
surv.formula <- Surv(survtime, status) ~ APOE4
tempfn <- function(response, RE = '(1 + time|id)', interval = c(-2, 2)){
  data <- newadni(response)
  long.formula <- as.formula(paste0(response, ' ~ ', paste('time','age_scaled','APOE4', sep = '+'), ' + ', RE))
  formulas <- parseFormula(long.formula)
  surv <- parseCoxph(surv.formula, data)
  n <- surv$n
  
  #' Initial conditions ----
  
  inits.long <- Longit.inits(long.formula, disp.formula, data)
  inits.surv <- TimeVarCox(data, inits.long$b, surv$ph, formulas)
  
  # Longitudinal parameters
  beta <- inits.long$beta.init
  D <- inits.long$D.init
  b <- lapply(1:n, function(i) inits.long$b[i, ])
  
  # Survival parameters
  zeta <- inits.surv$inits[match(colnames(surv$ph$x), names(inits.surv$inits))]
  names(zeta) <- paste0('zeta_', names(zeta))
  gamma <- inits.surv$inits[grepl('gamma', names(inits.surv$inits))]
  
  #' Data objects ----
  sv <- surv.mod(surv$ph, surv$survdata, formulas, inits.surv$l0.init)
  dmats <- createDataMatrices(data, formulas, disp.formula)
  # Truncation amount.
  summax <- max(sapply(dmats$Y, max)) * 2
  
  message('Doing optim...')
  optim.inits <- suppressMessages(get.delta.inits(dmats, beta, b, 'optim', summax, verbose = T, min.profile.length = 1,
                                                  interval = interval, percentile = c(.4,.6)))
  
  optim.inits$summax <- summax
  return(optim.inits)
}

# Command -----------------|  gen.summary | delta.optim |#   
tempfn('ADAS11')          #|        1.500 |   0.34/0.42 |#
tempfn('ADAS13')          #|        1.446 |   0.31/0.39 |#
# MMSE requires intercept-only RE structure; done below...
tempfn('RAVLT.immediate') #|        1.620 |   0.37/0.46 |#
tempfn('RAVLT.learning')  #|        1.222 |   0.35/0.26 |#
tempfn('RAVLT.forgetting')#|        2.111 |   0.88/0.64 |#
tempfn('FAQ')             #|        0.663 |  -0.31/0.45 |#

#' Candidate measures may then be ADAS(..), RAVLT for underdispersed and FAQ for overdispersed.
tempfn('MMSE', '(1|id)')  # Heavily underdispersed ...
                          # Might be worth expanding boundary conditions from [-2, 2] to [-6, 6] or something??
tempfn('MMSE', RE = '(1|id)', interval = c(-5, 5))
                          #|      11.742  |   2.43/2.50 |#






