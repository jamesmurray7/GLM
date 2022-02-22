data <- dd$data
newdata <- data %>% filter(id ==1, time < 4)
pu <- prepdata(newdata, id = 1, u = 4, fit = fit)
pt <- prepdata(newdata, id = 1, fit = fit)

pi.store <- numeric(100)
pb <- utils::txtProgressBar(max = 100, style = 3)

microbenchmark::microbenchmark(
  `b ~ N(b,S)` = {
    for(i in 1:100){
      O <- Omega.draw(fit)
      
      # mh <- b.mh(O, Sigmai.prop, b, pt)
      # b <- mh$b
      
      b <- b.draw(pt$b, # Both this and the MH way to draw from f(b) work, maybe see what's fastest?
                  pt$long$Xt, pt$long$Yt, pt$long$Zt,
                  O$beta, O$var.e, O$D,
                  pt$surv$Delta, pt$surv$K, pt$surv$Fi, pt$surv$l0i, pt$surv$KK.t,
                  pt$surv$Fu.t, pt$surv$l0u.t, O$gamma, O$eta)
      
      
      
      pi.store[i] <- S(b$b, O, pu$surv) / S(b$b, O, pt$surv)
      # utils::setTxtProgressBar(pb, i)
    }
  },
  `MH for f(b|...)` = {
    b <- pt$b
    Sigmai.prop <- pt$S
    for(i in 1:100){
      O <- Omega.draw(fit)
      
      mh <- b.mh(O, Sigmai.prop, b, pt)
      b <- mh$b
      
      pi.store[i] <- S(b, O, pu$surv) / S(b, O, pt$surv)
      utils::setTxtProgressBar(pb, i)
    }
    
    # Output median [95% CI] for survival at time u.
    quantile(pi.store, p = c(.025, .5, .975))
  },
  times = 25
) -> bench


plot(bench)
