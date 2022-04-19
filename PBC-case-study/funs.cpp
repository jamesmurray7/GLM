#include <RcppArmadillo.h>
#include <math.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


// Log-likelihoods --------------------------------------------------------
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

// 2. Scores for linear predictors ----------------------------------------
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


// 3. Defining a joint density --------------------------------------------
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
  }else if(family == "negative.binomial"){
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
  }else if(family == "negative.binomial"){
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

// 4. Defining updates to fixed-effects \beta -----------------------------

// [[Rcpp::export]]
vec Sbeta(vec& beta, mat& X, vec& Y, mat& Z, vec& b, double sigma, std::string& family){
  int p = beta.size();
  vec Score = vec(p);
  vec eta = X * beta + Z * b;
  
  if(family == "gaussian"){
    Score += X.t() * Score_eta_gauss(Y, eta, sigma);
  }else if(family == "binomial"){
    Score += X.t() * Score_eta_binom(Y, eta);
  }else if(family == "poisson"){
    Score += X.t() * Score_eta_poiss(Y, eta);
  }else if(family == "negative.binomial"){
    Score += X.t() * Score_eta_negbin(Y, eta, sigma);
  }
  
  return Score;
}

// [[Rcpp::export]]
mat Hbeta(vec& beta, mat& X, vec& Y, mat& Z, vec& b, double sigma, std::string& family, long double eps){
  int p = beta.size();
  mat out = zeros<mat>(p, p);
  vec S0 = Sbeta(beta, X, Y, Z, b, sigma, family);
  for(int i = 0; i < p; i++){
    vec bb = beta;
    double xi = std::max(beta[i], 1.0);
    bb[i] = beta[i] + eps * xi;
    vec Sdiff = Sbeta(bb, X, Y, Z, b, sigma, family) - S0;
    out.col(i) = Sdiff/(bb[i]-beta[i]);
  }
  return 0.5 * (out + out.t()); // Ensure symmetry
}

// 5. Defining updates to dispersion/scale parameters ---------------------

// 5a. Update for \theta, in case of Y|. ~ NB(mu, theta) parameterisation.
double ll_theta_quad(double theta, vec& beta, mat& X, vec& Y, mat& Z, vec& b, mat& S,
                     vec& w, vec& v){
  vec rhs = vec(Y.size());
  vec tau = sqrt(diagvec(Z * S * Z.t()));
  vec eta = X * beta + Z * b;
  for(int l = 0; l < w.size(); l++){
    vec this_eta = eta + v[l] * tau;
    rhs += w[l] * log(exp(this_eta) + theta);
  }
  vec out = vec(Y.size());
  for(int i = 0; i < Y.size(); i++){
    out[i] = lgamma(theta + Y[i]) - lgamma(theta) + theta * log(theta) - (theta + Y[i]) * rhs[i];
  }
  return sum(out);
}

// [[Rcpp::export]]
double Stheta(double theta, vec& beta, mat& X, vec& Y, mat& Z, vec& b, mat& S,
              vec& w, vec& v, double eps){
  double theta1 = theta + eps;
  double theta2 = theta - eps;
  double l1 = ll_theta_quad(theta1, beta, X, Y, Z, b, S, w, v);
  double l2 = ll_theta_quad(theta2, beta, X, Y, Z, b, S, w, v);
  return (l1-l2)/(theta1-theta2);
}

// [[Rcpp::export]]
double Htheta(double theta, vec& beta, mat& X, vec& Y, mat& Z, vec& b, mat& S,
              vec& w, vec& v, double eps){
  double f0 = Stheta(theta, beta, X, Y, Z, b, S, w, v, eps);
  double xi = eps * (abs(theta) + eps);
  double tt = theta + xi;
  double fdiff = Stheta(tt, beta, X, Y, Z, b, S, w, v, eps) - f0;
  return fdiff/(tt - theta);
}

// 5b. Update for residual variance
// [[Rcpp::export]]
long double vare_update(mat& X, vec& Y, mat& Z, vec& b, vec& beta, vec& tau, vec& w, vec& v){
  vec eta = X * beta + Z * b;
  int gh = w.size();
  double out = 0.0;
  for(int l = 0; l < gh; l++){
    vec rhs = Y - eta - tau * v[l];
    out += w[l] * as_scalar(rhs.t() * rhs);
  }
  return out;
}

// 6. Defining updates for the survival pair (gamma, zeta)
// Define the conditional expectation and then take Score AND Hessian via forward differencing
double Egammazeta(vec& gammazeta, vec& b, mat& Sigma,
                  rowvec& S, mat& SS, mat& Fu, rowvec& Fi, vec& haz, int Delta, vec& w, vec& v){
  double g = as_scalar(gammazeta.head(1)); // gamma will always be scalar and the first element
  vec z = gammazeta.subvec(1, gammazeta.size() - 1);  // with the rest of the vector constructed by zeta
  // determine tau
  vec tau = pow(g, 2.0) * diagvec(Fu * Sigma * Fu.t());
  double rhs = 0.0;
  for(int l = 0; l < w.size(); l++){
    rhs += w[l] * as_scalar(haz.t() * exp(SS * z + Fu * (g * b) + v[l] * pow(tau, 0.5)));
  }
  return as_scalar(Delta * (S * z + Fi * (g * b)) - rhs);
}
// [[Rcpp::export]]
vec Sgammazeta(vec& gammazeta, vec& b, mat& Sigma,
              rowvec& S, mat& SS, mat& Fu, rowvec& Fi, vec& haz, int Delta, vec& w, vec& v, long double eps){
  vec out = vec(gammazeta.size());
  double f0 = Egammazeta(gammazeta, b, Sigma, S, SS, Fu, Fi, haz, Delta, w, v);
  for(int i = 0; i < gammazeta.size(); i++){
    vec ge = gammazeta;
    double xi = std::max(ge[i], 1.0);
    ge[i] = gammazeta[i] + xi * eps;
    double fdiff = Egammazeta(ge, b, Sigma, S, SS, Fu, Fi, haz, Delta, w, v) - f0;
    out[i] = fdiff/(ge[i]-gammazeta[i]);
  }
  return out;
}
// [[Rcpp::export]]
mat Hgammazeta(vec& gammazeta, vec& b, mat& Sigma,
              rowvec& S, mat& SS, mat& Fu, rowvec& Fi, vec& haz, int Delta, vec& w, vec& v, long double eps){
  mat out = zeros<mat>(gammazeta.size(), gammazeta.size());
  vec f0 = Sgammazeta(gammazeta, b, Sigma, S, SS, Fu, Fi, haz, Delta, w, v, eps);
  for(int i = 0; i < gammazeta.size(); i++){
    vec ge = gammazeta;
    double xi = std::max(ge[i], 1.0);
    ge[i] = gammazeta[i] + xi * eps;
    vec fdiff = Sgammazeta(ge, b, Sigma, S, SS, Fu, Fi, haz, Delta, w, v, eps) - f0;
    out.col(i) = fdiff/(ge[i]-gammazeta[i]);
  }
  return 0.5 * (out + out.t());
}

// 7. Defining update to the baseline hazard lambda_0.
// [[Rcpp::export]]
mat lambdaUpdate(List survtimes, mat& ft, double gamma, vec& zeta,
                 List S, List Sigma, List b, vec& w, vec& v){
  int gh = w.size();
  int id = b.size();
  mat store = zeros<mat>(ft.n_rows, id);
  for(int i = 0; i < id; i++){
    vec survtimes_i = survtimes[i];
    mat Sigma_i = Sigma[i];
    vec b_i = b[i];
    rowvec S_i = S[i];
    for(int j = 0; j < survtimes_i.size(); j++){
      rowvec Fst = ft.row(j);
      double tau = as_scalar(pow(gamma, 2.0) * Fst * Sigma_i * Fst.t());
      vec rhs = gamma * b_i; //vec(b_i.size());
      double mu = as_scalar(exp(S_i * zeta + Fst * rhs));
      for(int l = 0; l < gh; l++){
        store(j, i) += as_scalar(w[l] * mu * exp(v[l] * sqrt(tau)));
      }
    }
  }
  return store;
}

