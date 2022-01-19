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
vec calcZ(vec& lambda, vec& nu, int summax){
  vec out = vec(lambda.size());
  vec js = SEQ_Z(summax);
  for(int i = 0; i < out.size(); i++){
    for(int j = 0; j < js.size(); j++){
      out[i] += pow(lambda[i], js[j])/pow(gamma(js[j] + 1.0), nu[i]);
    }
  }
  return out;
}

// [[Rcpp::export]]
double calcZ_scalar(double lambda, double nu, int summax){
  double out = 0.0;
  vec js = SEQ_Z(summax);
  for(int j = 0; j < js.size(); j++){
    out += pow(lambda, js[j])/pow(gamma(js[j] + 1.0), nu);
  }
  return out;
}

// [[Rcpp::export]]
double mu_lambda_eq(double lambda, double mu, double nu, int summax){
  vec js = SEQ_Z(summax);
  double out = 0.0;
  for(int j = 0; j < js.size(); j++){
    out += ((js[j] - mu) * pow(lambda, js[j]))/pow(gamma(js[j] + 1.0), nu);
  }
  return out;
}

// [[Rcpp::export]]
double mu_lambdaZ_eq(double lambda, double mu, double nu, int summax){
  vec js = SEQ_Z(summax);
  double Z = calcZ_scalar(lambda, nu, summax);
  double rhs = 0.0;
  for(int j = 0; j < js.size(); j++){
    rhs += (js[j] * pow(lambda, js[j])) / (pow(gamma(js[j] + 1.0), nu) * Z);
  }
  return mu - rhs;
}

// [[Rcpp::export]]
vec cmp_pmf(vec& Y, vec& lambda, vec& nu, int summax){
  vec out = vec(lambda.size());
  vec Z = calcZ(lambda, nu, summax);
  vec l = Y % log(lambda) - nu % lgamma(Y + 1.0) - log(Z);
  return exp(l);
}

// [[Rcpp::export]]
vec cmp_pmf_scalar(vec& Y, double lambda, double nu, int summax){
  vec out = vec(Y.size());
  double Z = calcZ_scalar(lambda, nu, summax);
  vec L = exp(Y * log(lambda) - nu * lgamma(Y + 1.0) - log(Z));
  L.replace(datum::inf, 1e100);
  return L;
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

