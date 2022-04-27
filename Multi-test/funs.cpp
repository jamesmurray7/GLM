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
static double log2pi = log(2.0 * M_PI); // nice to know this works?

// [[Rcpp::export]]
double joint_density(vec& b, List Y, List X, List Z,                  // Longitudinal + Random effects.
                     vec& beta, mat& D, List sigma, List family,
                     int Delta, rowvec& S, rowvec& Fi, double l0i,
                     mat& SS, mat& Fu, rowvec& haz, vec& gamma_rep, vec& zeta,
                     List beta_inds, List b_inds, int K){
  double ll = 0.0; // Compile longitudinal parts ---------
  for(int k = 0; k < K; k++){
    vec Yk = Y[k];
    mat Xk = X[k];
    mat Zk = Z[k];
    std::string f = family[k];
    double sigmak = sigma[k];
    uvec beta_k_inds = beta_inds[k];
    uvec b_k_inds = b_inds[k];
    vec beta_k = beta.elem(beta_k_inds); // Ensure indexing from zero!!
    vec b_k = b.elem(b_k_inds);
    vec eta = Xk * beta_k + Zk * b_k;
    // Rcout << "k: " << k << std::endl;
    // Rcout << "Family k: " << f << std::endl;
    // Rcout << "eta_k: " << eta << std::endl;
    if(f == "gaussian"){
      ll += gaussian_ll(Yk, eta, sigmak);
    }else if(f == "binomial"){
      ll += binomial_ll(Yk, eta);
    }else if(f == "poisson"){
      ll += poisson_ll(Yk, eta);
    }else if(f == "negative.binomial"){
      ll += negbin_ll(Yk, eta, sigmak);
    }
  }
  // Rcout << "ll: " << ll << std::endl;
  int q = b.size();
  double temp = 0.0;
  if(Delta == 1) temp = log(l0i);
  return -ll + -1.0 * as_scalar(-q/2.0 * log2pi - 0.5 * log(det(D)) - 0.5 * b.t() * D.i() * b + 
                                temp + Delta * (S * zeta + Fi * (gamma_rep % b)) - haz * exp(SS * zeta + Fu * (gamma_rep % b)));
}

// First derivative of the joint density with respect to b.
// [[Rcpp::export]]
vec joint_density_ddb(vec& b, List Y, List X, List Z,                  // Longitudinal + Random effects.
                      vec& beta, mat& D, List sigma, List family,
                      int Delta, rowvec& S, rowvec& Fi, double l0i,
                      mat& SS, mat& Fu, rowvec& haz, vec& gamma_rep, vec& zeta,
                      List beta_inds, List b_inds, int K){
  vec Score = vec(b.size());
  for(int k = 0; k < K; k++){
    vec Yk = Y[k];
    mat Xk = X[k];
    mat Zk = Z[k];
    std::string f = family[k];
    double sigmak = sigma[k];
    uvec beta_k_inds = beta_inds[k];
    uvec b_k_inds = b_inds[k];
    vec beta_k = beta.elem(beta_k_inds); // Ensure indexing from zero!!
    vec b_k = b.elem(b_k_inds);
    vec eta = Xk * beta_k + Zk * b_k;
    int qk = b.size();
    if(f == "gaussian"){
      Score.elem(b_k_inds) += Zk.t() * Score_eta_gauss(Yk, eta, sigmak);
    }else if(f == "binomial"){
      Score.elem(b_k_inds) += Zk.t() * Score_eta_binom(Yk, eta);
    }else if(f == "poisson"){
      Score.elem(b_k_inds) += Zk.t() * Score_eta_poiss(Yk, eta);
    }else if(f == "negative.binomial"){
      Score.elem(b_k_inds) += Zk.t() * Score_eta_negbin(Yk, eta, sigmak);
    }
  }
  return -Score + -1.0 * (-D.i() * b + Delta * (Fi.t() % gamma_rep) - 
                               gamma_rep % (Fu.t() * (haz.t() % exp(SS * zeta + Fu * (gamma_rep % b)))));
}

// Second derivative of the joint density with respect to b. This is obtained via forward differencing.
// NB sometimes this doesn't work, and hessian obtained from optim (BFGS) is taken instead (as yet unknown why they're different).
// [[Rcpp::export]]
mat joint_density_sdb(vec& b, List Y, List X, List Z,                  // Longitudinal + Random effects.
                      vec& beta, mat& D, List sigma, List family,
                      int Delta, rowvec& S, rowvec& Fi, double l0i,
                      mat& SS, mat& Fu, rowvec& haz, vec& gamma_rep, vec& zeta,
                      List beta_inds, List b_inds, int K, double eps){
  int n = b.size();
  mat out = zeros<mat>(n, n);
  vec f0 = joint_density_ddb(b, Y, X, Z, beta, D, sigma, family, Delta, S, Fi, l0i, SS, Fu, haz, gamma_rep, zeta, beta_inds, b_inds, K);
  for(int i = 0; i < n; i++){
    vec bb = b;
    double xi = std::max(b[i], 1.0);
    bb[i] = b[i] + (eps * xi);
    vec fdiff = joint_density_ddb(bb, Y, X, Z, beta, D, sigma, family, Delta, S, Fi, l0i, SS, Fu, haz, gamma_rep, zeta, beta_inds, b_inds, K);
    out.col(i) = fdiff/(bb[i]-b[i]);
  }
  return 0.5 * (out + out.t()); // Ensure symmetry
}

// 4. Defining updates to fixed-effects \beta -----------------------------

// [[Rcpp::export]]
vec Sbeta(vec& beta, List X, List Y, List Z, List b, List sigma, List family, List beta_inds, int K){
  int p = beta.size();
  vec Score = vec(p);
  
  for(int k = 0; k < K; k++){
    vec Yk = Y[k];
    mat Xk = X[k];
    mat Zk = Z[k];
    std::string f = family[k];
    double sigmak = sigma[k];
    uvec beta_k_inds = beta_inds[k];
    vec beta_k = beta.elem(beta_k_inds); // Ensure indexing from zero!!
    vec b_k = b[k];
    vec eta = Xk * beta_k + Zk * b_k;
    int qk = b.size();
    if(f == "gaussian"){
      Score.elem(beta_k_inds) += Xk.t() * Score_eta_gauss(Yk, eta, sigmak);
    }else if(f == "binomial"){
      Score.elem(beta_k_inds) += Xk.t() * Score_eta_binom(Yk, eta);
    }else if(f == "poisson"){
      Score.elem(beta_k_inds) += Xk.t() * Score_eta_poiss(Yk, eta);
    }else if(f == "negative.binomial"){
      Score.elem(beta_k_inds) += Xk.t() * Score_eta_negbin(Yk, eta, sigmak);
    }
  }
  
  return Score;
}

// [[Rcpp::export]]
mat Hbeta(vec& beta, List X, List Y, List Z, List b, List sigma, List family, List beta_inds, int K, double eps){
  int p = beta.size();
  mat out = zeros<mat>(p, p);
  vec S0 = Sbeta(beta, X, Y, Z, b, sigma, family, beta_inds, K);
  for(int i = 0; i < p; i++){
    vec bb = beta;
    double xi = std::max(beta[i], 1.0);
    bb[i] = beta[i] + eps * xi;
    vec Sdiff = Sbeta(bb, X, Y, Z, b, sigma, family, beta_inds, K) - S0;
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
double Egammazeta(vec& gammazeta, vec& b, List Sigma,
                  rowvec& S, mat& SS, mat& Fu, rowvec& Fi, vec& haz, int Delta, vec& w, vec& v,
                  List b_inds, int K, int q){
  vec g = gammazeta.head(K); // First K elements as proportional association only.
  vec z = gammazeta.subvec(K, gammazeta.size() - 1);  // with the rest of the vector constructed by zeta
  // determine tau
  vec tau = vec(Fu.n_rows);
  vec gammas = vec(q);
  for(int k = 0; k < K; k++){
    double gk = g[k];
    uvec b_inds_k = b_inds[k];
    gammas.elem(b_inds_k) += gk;
    mat Fu_k = Fu.cols(b_inds_k);
    mat Sigma_k = Sigma[k];
    tau += pow(gk, 2.0) * diagvec(Fu_k * Sigma_k * Fu_k.t());
  }
  double rhs = 0.0;
  for(int l = 0; l < w.size(); l++){
    rhs += w[l] * as_scalar(haz.t() * exp(SS * z + Fu * (b % gammas) + v[l] * pow(tau, 0.5)));
  }
  return as_scalar(Delta * (S * z + Fi * (b % gammas)) - rhs);
}
// [[Rcpp::export]]
vec Sgammazeta(vec& gammazeta, vec& b, List Sigma,
               rowvec& S, mat& SS, mat& Fu, rowvec& Fi, vec& haz, int Delta, vec& w, vec& v,
               List b_inds, int K, int q, long double eps){
  vec out = vec(gammazeta.size());
  double f0 = Egammazeta(gammazeta, b, Sigma, S, SS, Fu, Fi, haz, Delta, w, v, b_inds, K, q);
  for(int i = 0; i < gammazeta.size(); i++){
    vec ge = gammazeta;
    double xi = std::max(ge[i], 1.0);
    ge[i] = gammazeta[i] + xi * eps;
    double fdiff = Egammazeta(ge, b, Sigma, S, SS, Fu, Fi, haz, Delta, w, v, b_inds, K, q) - f0;
    out[i] = fdiff/(ge[i]-gammazeta[i]);
  }
  return out;
}
// [[Rcpp::export]]
mat Hgammazeta(vec& gammazeta, vec& b, List Sigma,
               rowvec& S, mat& SS, mat& Fu, rowvec& Fi, vec& haz, int Delta, vec& w, vec& v,
               List b_inds, int K, int q, long double eps){
  mat out = zeros<mat>(gammazeta.size(), gammazeta.size());
  vec f0 = Sgammazeta(gammazeta, b, Sigma, S, SS, Fu, Fi, haz, Delta, w, v, b_inds, K, q, eps);
  for(int i = 0; i < gammazeta.size(); i++){
    vec ge = gammazeta;
    double xi = std::max(ge[i], 1.0);
    ge[i] = gammazeta[i] + xi * eps;
    vec fdiff = Sgammazeta(ge, b, Sigma, S, SS, Fu, Fi, haz, Delta, w, v, b_inds, K, q, eps) - f0;
    out.col(i) = fdiff/(ge[i]-gammazeta[i]);
  }
  return 0.5 * (out + out.t());
}

// 7. Defining update to the baseline hazard lambda_0.
// [[Rcpp::export]]
mat lambdaUpdate(List survtimes, mat& ft, vec& gamma, vec& gamma_rep, vec& zeta,
                 List S, List Sigma, List b, vec& w, vec& v, List b_inds, int K, int q){
  int gh = w.size();
  int n = b.size();
  mat store = zeros<mat>(ft.n_rows, n);
  for(int i = 0; i < n; i++){ // Loop over i subjects
    vec survtimes_i = survtimes[i];
    List Sigma_i = Sigma[i];
    vec b_i = b[i];
    rowvec S_i = S[i];
    for(int j = 0; j < survtimes_i.size(); j++){ // Loop over subject i's j survived failure times.
      rowvec Fst = ft.row(j);
      double tau = 0.0;
      vec rhs = gamma_rep % b_i;
      for(int k = 0; k < K; k++){ // Loop over the K longitudinal responses.
        mat Sigma_ik = Sigma_i[k];
        uvec b_inds_k = b_inds[k];
        rowvec Fst_k = Fst.elem(b_inds_k).t();
        tau += as_scalar(pow(gamma[k], 2.0) * Fst_k * Sigma_ik * Fst_k.t());
      }
      double mu = as_scalar(exp(S_i * zeta + Fst * rhs));
      for(int l = 0; l < gh; l++){
        store(j, i) += as_scalar(w[l] * mu * exp(v[l] * sqrt(tau)));
      }
    }
  }
  return store;
}
