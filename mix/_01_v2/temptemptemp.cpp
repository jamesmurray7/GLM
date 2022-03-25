#include <RcppArmadillo.h>
#include <math.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
vec lfact(vec& v){
  vec out = vec(v.size());
  for(int i = 0; i < v.size(); i++){
    out[i] = lgamma(v[i] + 1.0);
  }
  return out;
}

// Log-likelihood for Poisson
// [[Rcpp::export]]
double poisson_ll(vec& Y, vec& eta){
  vec out = vec(Y.size());
  for(int i = 0; i < Y.size(); i++){
    out[i] = R::dpois(Y[i], exp(eta[i]), 1);
  }
  return sum(out);
}

// With quadrature

// [[Rcpp::export]]
double poisson_ll_quad(vec& beta, mat& X, vec& Y, mat& Z, vec& b,
                       mat& S, vec& w, vec& v){
  vec eta = X * beta + Z * b;
  Rcout << "eta" << eta << std::endl;
  vec lfactY = lfact(Y);
  vec out = vec(Y.size()), rhs = out;
  // quadrature
  int gh = w.size();
  vec tau = sqrt(diagvec(Z * S * Z.t()));
  Rcout << "tau" << tau << std::endl;
  for(int l = 0; l < gh; l++){
    vec this_eta = eta + v[l] * tau;
    rhs += w[l] * exp(this_eta);
  }
  out = Y % eta - rhs - lfactY;
  return sum(out);
}

