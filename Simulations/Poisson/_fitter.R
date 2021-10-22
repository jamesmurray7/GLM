files <- dir(pattern = 'sim\\d.RData', path="~/Downloads")

files <- as.list(files)

fn <- function(f){
  message(f)
  load(paste0("~/Downloads/", f))
  pb <- utils::txtProgressBar(max=50, style = 3)
  fits <- list()
  for(m in 1:50){
    dat <- full.data[[m]]
    ph <- coxph(Surv(survtime, status) ~ cont + bin, data = dplyr::distinct(dat, id, cont, bin, survtime, status))
    fits[[m]] <- tryCatch(
      suppressMessages(em(dat, ph, gh.nodes = 3, nK = 2)),
      error = function(e) NULL
    )
    setTxtProgressBar(pb, m)
  }
  out.destination <- paste0(getwd())
  message("\nSaving fits for", f, "in ", out.destination, 'fits3-', f)
  save(fits, file = paste0(out.destination, "fits3-", f))
}

lapply(files, fn)

fn <- function(f){
  message(f)
  load(paste0("~/Downloads/", f))
  pb <- utils::txtProgressBar(max=50, style = 3)
  fits <- list()
  for(m in 1:50){
    dat <- full.data[[m]]
    ph <- coxph(Surv(survtime, status) ~ cont + bin, data = dplyr::distinct(dat, id, cont, bin, survtime, status))
    fits[[m]] <- tryCatch(
      suppressMessages(em(dat, ph, gh.nodes = 9, nK = 2)),
      error = function(e) NULL
    )
    setTxtProgressBar(pb, m)
  }
  out.destination <- paste0(getwd())
  message("\nSaving fits for", f, "in ", out.destination, 'fits9-', f)
  save(fits, file = paste0(out.destination, "fits9-", f))
}

lapply(files, fn)
