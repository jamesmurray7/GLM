// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <math.h>
// #include <float.h>
using namespace Rcpp;
using namespace arma;

vec psigamma(vec & x, int deriv){
  int n = x.size();
  vec xx = vec(n);
  for(int a = 0; a < n ; a++){
    xx[a] = R::psigamma(x[a], deriv);
  }
  return xx;
}

// [[Rcpp::export]]
double dgamma_(const vec& Y, 
              const double& shape, const vec& scale){

  /* ****
   * Shape: alpha (double, phi);
   * Scale: beta  (vector, mu/phi)
   * **** */
  vec frac = (1./(pow(scale, shape) * tgamma(shape)));
  vec out = frac % pow(Y, shape - 1.) % exp(-Y/scale);
  return sum(log(out));
}

// [[Rcpp::export]]
