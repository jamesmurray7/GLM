#include <RcppArmadillo.h>
#include <math.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


// Log-likelihoods -------------------------------------
// Gaussian
// [[Rcpp::export]]
double gaussian_ll(vec& Y, vec& eta, double sigma){ // important/misleading - sigma is the __variance__!!
  vec out = vec(Y.size());
  for(int i = 0; i < Y.size(); i++){
    out[i] = R::dnorm(Y[i], eta[i], sqrt(sigma), 1);
  }
  return sum(out);
}

// Binomial
// [[Rcpp::export]]
double binomial_ll(vec& Y, vec& eta){
  vec mu = exp(eta)/(1.0 + exp(eta));
  vec out = vec(mu.size());
  for(int i = 0; i < mu.size(); i++){
    out[i] = R::dbinom(Y[i], 1, mu[i], 1);
  }
  return sum(out);
}

// Poisson
// [[Rcpp::export]]
double poisson_ll(vec& Y, vec& eta){
  vec out = vec(Y.size());
  for(int i = 0; i < Y.size(); i++){
    out[i] = R::dpois(Y[i], exp(eta[i]), 1);
  }
  return sum(out);
}

// Negative Binomial
// [[Rcpp::export]]
double negbin_ll(vec& Y, vec& eta, double theta){
  vec out = vec(Y.size());
  for(int i = 0; i < Y.size(); i++){
    out[i] = R::dnbinom_mu(Y[i], theta, exp(eta[i]), 1);
  }
  return sum(out);
}

// 2. Scores for linear predictors -------------------------------------
// Gaussian -> binom -> poisson -> negbin
// [[Rcpp::export]]
vec Score_eta_gauss(vec& Y, vec& eta, double sigma){
  int m = Y.size();
  mat V = mat(m, m, fill::eye);
  V *= sigma;
  return V.i() * (Y - eta);
}
// [[Rcpp::export]]
vec Score_eta_binom(vec& Y, vec& eta){
  return Y - exp(eta) / (exp(eta) + 1.0);
}
// [[Rcpp::export]]
vec Score_eta_poiss(vec& Y, vec& eta){
  return Y - exp(eta);
}
// [[Rcpp::export]]
vec Score_eta_negbin(vec& Y, vec& eta, double theta){
  return (theta * (Y - exp(eta))) / (exp(eta) + theta);
}


// 3. Defining a joint density ------------------------------------------
static double log2pi = log(2.0 * M_PI); // does this work?

// [[Rcpp::export]]
double joint_density(vec& b, vec& Y, mat& X, mat& Z,      // Longitudinal + REs
                     vec& beta, mat& D, double sigma, std::string& family,
                     int Delta, rowvec& S, rowvec& Fi, double l0i,  // Everything else.
                     mat& SS, mat& Fu, rowvec& haz, double gamma, vec& zeta){
  vec eta = X * beta + Z * b;
  
  // Determine longitudinal response log-likelihood.
  double ll = 0.0;
  if(family == "gaussian"){
    ll += gaussian_ll(Y, eta, sigma);
  }else if(family == "binomial"){
    ll += binomial_ll(Y, eta);
  }else if(family == "poisson"){
    ll += poisson_ll(Y, eta);
  }else if(family == "negative binomial"){
    ll += negbin_ll(Y, eta, sigma);
  }
  
  int q = b.size();
  double temp = 0.0;
  if(Delta == 1) temp = log(l0i);
  return -ll + -1.0 * as_scalar(-q/2.0 * log2pi - 0.5 * log(det(D)) - 0.5 * b.t() * D.i() * b + 
                                  temp + Delta * (S * zeta + Fi * (gamma * b)) - haz * exp(SS * zeta + Fu * (gamma * b)));
}

// First derivative of the joint density with respect to b.
// [[Rcpp::export]]
vec joint_density_ddb(vec& b, vec& Y, mat& X, mat& Z,      // Longitudinal + REs
                      vec& beta, mat& D, double sigma, std::string& family,
                      int Delta, rowvec& S, rowvec& Fi, double l0i,  // Everything else.
                      mat& SS, mat& Fu, rowvec& haz, double gamma, vec& zeta){
  int q = b.size();
  vec eta = X * beta + Z * b;
  
  // Determine longitudinal response log-likelihood.
  vec Score = vec(q);
  if(family == "gaussian"){
    Score += Z.t() * Score_eta_gauss(Y, eta, sigma);
  }else if(family == "binomial"){
    Score += Z.t() * Score_eta_binom(Y, eta);
  }else if(family == "poisson"){
    Score += Z.t() * Score_eta_poiss(Y, eta);
  }else if(family == "negative binomial"){
    Score += Z.t() * Score_eta_negbin(Y, eta, sigma);
  }
  
  return -Score + -1.0 * (-D.i() * b + Delta * gamma * Fi.t() - 
                               gamma * Fu.t() * (haz.t() % exp(SS * zeta + gamma * Fu * b)));
}

// Second derivative of the joint density with respect to b. This is obtained via forward differencing.
// [[Rcpp::export]]
mat joint_density_sdb(vec& b, vec& Y, mat& X, mat& Z,      // Longitudinal + REs
                      vec& beta, mat& D, double sigma, std::string& family,
                      int Delta, rowvec& S, rowvec& Fi, double l0i,  // Everything else.
                      mat& SS, mat& Fu, rowvec& haz, double gamma, vec& zeta, long double eps){
  int n = b.size();
  mat out = zeros<mat>(n, n);
  vec f0 = joint_density_ddb(b, Y, X, Z, beta, D, sigma, family, Delta, S, Fi, l0i, SS, Fu, haz, gamma, zeta);
  for(int i = 0; i < n; i++){
    vec bb = b;
    double xi = std::max(b[i], 1.0);
    bb[i] = b[i] + (eps * xi);
    vec fdiff = joint_density_ddb(bb, Y, X, Z, beta, D, sigma, family, Delta, S, Fi, l0i, SS, Fu, haz, gamma, zeta);
    out.col(i) = fdiff/(bb[i]-b[i]);
  }
  return 0.5 * (out + out.t()); // Ensure symmetry
}

