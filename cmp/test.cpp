#include <RcppArmadillo.h>
#include <math.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::vec SEQ(long double from, long double to, long unsigned int length_out){
  return arma::linspace(from, to, length_out);
}

// [[Rcpp::export]]
arma::vec SEQ_Z(long double summax){ // Shorthand for Z_(lambda_i, nu_i)
  return arma::linspace(0, summax, summax + 1.0);
}

// [[Rcpp::export]]
arma::vec SEQ_Z_2(long double summax){ // Shorthand for Z_(lambda_i, nu_i)
  return arma::linspace(1, summax, summax);
}

// [[Rcpp::export]]
vec calc_logZ(vec& eta, vec& nu, int summax){
  vec ns  = SEQ_Z(summax); // Initialise vector of ns (n = 0, ..., summax)
  vec out = vec(eta.size());
  for(int i = 0; i < eta.size(); i ++){
    for(int j = 1; j < ns.size(); j++){
      out[i] += ns[j] * eta[i] - lgamma(ns[j] + 1.0) * nu[i];
    }
    
  }
  return out + 1.0; // from cumsum log(0)
}  

// [[Rcpp::export]]
vec calc_logZ2(vec& eta, vec& nu, int summax){
  vec ns  = SEQ_Z_2(summax); // Initialise vector of ns (n = 1, ..., summax)
  vec lcs = cumsum(log(ns)); // log cumsum
  vec out = vec(eta.size());
  for(int i = 0; i < eta.size(); i ++){
    out[i] = 1.0 + sum(ns * eta[i] - lcs * nu[i]);
  }
  return out; // from cumsum log(0)
}  

// vec dCMP(vec& Y, )

// rcomp from mpcmp package
//for (i in 1:n) {
    //if (mu[i] == 0 | lambda[i] == 0) {
      //x[i] <- 0
    //} else if (mu[i] < 0 | lambda[i] < 0 | nu[i] <= 0) {
      //x[i] <- NA
      //warn <- TRUE
    //} else {
      //y <- 0
      //dc <- dcomp(0:summax, nu = nu[i], lambda = lambda[i], summax = summax)
      //py <- dc[y + 1]
      //while (py <= unif[i]) {
        //y <- y + 1
        //py <- py + dc[y + 1]
      //}
      //x[i] <- y
    //}
  //}

// dcomp from mpcmp package...
  //pmf <- rep(0, length(x))
  //for (i in 1:length(x)) {
    //if ((mu[i] == 0 || lambda[i] == 0) && x[i] == 0) {
      //pmf[i] <- 0 # log(1), 1 as the distribution is degenerated at 0
    //} else if (mu[i] < 0 | lambda[i] < 0 | nu[i] <= 0) {
      //pmf[i] <- NaN
      //warn <- TRUE
    //} else {
      //if (!is.wholenumber(x[i])) {
        //warning(paste("non-integer x =", x[i]))
        //pmf[i] <- -Inf # log(0)
      //} else {
        //if (x[i] < 0) {
          //pmf[i] <- -Inf
        //} else { # log(0)
          //# pmf <- log(density)
          //pmf[i] <- x[i] * log(lambda[i]) - (nu[i] * lfactorial(x[i])) -
            //logZ(log(lambda[i]), nu[i], summax)
        //}
      //}
    //}
  //}
  //if (!log.p) {
    //pmf <- exp(pmf)
  //}
  //return(pmf)
