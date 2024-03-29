files <- dir(pattern = 'sim', path="~/Downloads")
setwd('./Poisson/')
source('EM.R')
sourceCpp('../beta_test.cpp')
sourceCpp('../temp-gammaCalc.cpp')

files <- as.list(files)

fn <- function(f){
  message(f)
  load(paste0("~/Downloads/", f))
  pb <- utils::txtProgressBar(max=50, style = 3)
  fits3 <- fits9 <- list()
  for(m in 1:50){
    dat <- data[[m]]
    ph <- coxph(Surv(survtime, status) ~ cont + bin, data = dplyr::distinct(dat, id, cont, bin, survtime, status))
    fits3[[m]] <- tryCatch(
      suppressMessages(em(dat, ph, gh.nodes = 3, nK = 2)),
      error = function(e) NULL
    )
    fits9[[m]] <- tryCatch(
      suppressMessages(em(dat, ph, gh.nodes = 9, nK = 2)),
      error = function(e) NULL
    )
    setTxtProgressBar(pb, m)
  }
  out.destination <- paste0(getwd())
  message("\nSaving fits for", f, "in ", out.destination, '/fits3-', f)
  save(fits3, file = paste0(out.destination, "/fits3-", f))
  message("\nSaving fits for", f, "in ", out.destination, '/fits9-', f)
  save(fits9, file = paste0(out.destination, "/fits9-", f))
}

lapply(files, fn)

