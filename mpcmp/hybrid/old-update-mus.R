# Add to CURRENT dispersion estimates given potential new mu values.
all.mus.old <- all.mus
mm <- do.call(c, update$mus); 
.min.new <- max(0, min(mm) - 5); .max.new <- max(mm) + 5
all.mus <- generate_mus(.min.new,  .max.new);
.diff <- setdiff(all.mus, all.mus.old)

if(length(.diff) > 0){
  ld <- length(.diff)
  if(verbose) cat(sprintf("%d new mu elements.\n" , ld))
  
  # New lambda values for the new mu values...
  newmus <- generate_mus(.diff[1], .diff[ld])
  new.lambdas <- structure(gen_lambda_mat(mu_lower = .diff[1], mu_upper = .diff[ld],
                                          nus = nu.vec, summax = summax),
                           dimnames = list(as.character(newmus), as.character(nu.vec)))
  
  # New logZ values
  new.logZ <- structure(gen_logZ_mat(mu_lower = .diff[1], mu_upper = .diff[ld],
                                     nus = nu.vec, lambdamat = new.lambdas, summax = summax),
                        dimnames = list(as.character(newmus), as.character(nu.vec)))
  # New V values
  new.V <- structure(gen_V_mat(mu_lower = .diff[1], mu_upper = .diff[ld],
                               nus = nu.vec, lambdamat = new.lambdas, logZmat = new.logZ, summax = summax),
                     dimnames = list(as.character(newmus), as.character(nu.vec)))
  
  # Update lookup matrices mats for new values of mu...
  if(all(.diff > max(all.mus.old))){        # Case 1: If all new mu values are greater, then simply rbind below...
    if (verbose) cat('All new elements > old mus.\n')
    lambda.mat.new <- rbind(lambda.mat, new.lambdas)
    logZ.mat.new <- rbind(logZ.mat, new.logZ)
    V.mat.new <- rbind(V.mat, new.V)
  }else if(all(.diff < min(all.mus.old))){  # Case 2: If all new mu values are smaller, simply rbind above...
    if (verbose) cat('All new elements < old mus.\n')
    lambda.mat.new <- rbind(new.lambdas, lambda.mat)
    logZ.mat.new <- rbind(new.logZ, logZ.mat)
    V.mat.new <- rbinD(new.V, V.mat)
  }else{                             # Case 3: Mixture above and below old mus.
    lts <- which(.diffs < min(all.mus.old)); gts <- which(.diffs > max(all.mus.old))
    if (verbose) cat(sprintf("%d elements < old mus; %d elements > old mus.\n", lts, gts))
    lambda.mat.new <- rbind(rbind(new.lambdas[lts,,drop=F], lambda.mat), new.lambas[gts,,drop=F])
    logZ.mat.new <- rbind(rbind(new.logZ[lts,,drop=F], logZ.mat), new.logZ[gts,,drop=F])
    V.mat.new <- rbind(rbind(new.V[lts,,drop=F], V.mat), new.V[gts,,drop=F])
  }
}else{
  if(verbose) cat('No new mu elements to estimate.\n')
  lambda.mat.new <- lambda.mat; logZ.mat.new <- logZ.mat; V.mat.new <- V.mat
  all.mus <- all.mus.old # An exact subset, so just keep old values
  .min.new <- .min; .max.new <- .max
}