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
double ll_Gamma(const vec& Y, const double& shape, const vec& scale){
  int n = Y.size();
  vec out(n);
  for(int i = 0; i < n; i++){
    out[i]  = R::dgamma(Y[i], shape, scale[i], 1);
  }
  return sum(out);
}

// [[Rcpp::export]]
double ll_Gamma2(const vec& Y, const double& shape, const vec& mu){
  int n = Y.size();
  vec out = (shape - 1.) * log(Y) - lgamma(shape) - shape * log(mu) + 
                shape * log(shape) - (shape * Y)/mu;
  return sum(out);
}

// vec frac = (1./(pow(scale, shape) * tgamma(shape)));
// vec out = frac % pow(Y, shape - 1.) % exp(-Y/scale);
// return sum(log(out))


// log f(T_i, \Delta_i|\b; \Omega).
// [[Rcpp::export]]
double logfti(vec& b, rowvec& S, mat& SS, rowvec& Fi, mat& Fu,
              double l0i, rowvec& haz, int Delta, double gamma, vec& zeta){
  double temp = 0.0;
  if(Delta == 1) temp = log(l0i);
  
  return as_scalar(
    temp + Delta * (S * zeta + Fi * (gamma * b)) -
      haz * exp(SS * zeta + Fu * (gamma * b))
  );
}

// Defining a joint density
static double log2pi = log(2.0 * M_PI);
// The full joint density f(b,Y,T,...) and its gradient wrt. b
// [[Rcpp::export]]
double joint_density(vec& b, mat& X, vec& Y, mat& Z, vec& beta, double shape, 
                     mat& D, rowvec& S, mat& SS, rowvec& Fi, mat& Fu, double l0i, rowvec& haz,
                     int Delta, double gamma, vec& zeta){
  // Define mu, nu
  vec mu = exp(X * beta + Z * b);
  // Constituent parts.
  double ll_G = ll_Gamma(Y, shape, mu/shape);
  double ll_surv = logfti(b, S, SS, Fi, Fu, l0i, haz, Delta, gamma, zeta);
  int q = b.size();
  return -1. * (ll_G + as_scalar(-(double)q/2.0 * log2pi - 0.5 * log(det(D)) - 0.5 * b.t() * D.i() * b) + ll_surv);
}


// Derivatives
vec get_long_score(const vec& b, const mat& X, const vec& Y, const mat& Z, const vec& beta,
                   const double& shape, const mat& design){
  int q = design.n_cols;
  vec mu = exp(X * beta + Z * b), grad = vec(q);
  for(int qq = 0; qq < q; qq++){
    vec x = design.col(qq);
    grad[qq] = sum(shape * x % Y/mu - shape * x);
  }
  return grad;
}

mat get_long_hess(const vec& b, const mat& X, const vec& Y, const mat& Z, const vec& beta,
                  const double& shape, const mat& design){
  int mi = design.n_rows, q = design.n_cols;
  vec mu = exp(X * beta + Z * b);
  mat H = zeros<mat>(q, q);
  for(int j = 0; j < mi; j++){
    rowvec xjT = design.row(j);
    vec xj = xjT.t();
    H += (-shape * Y[j] / mu[j]) * xj * xjT;
  }
  return H;
}

// [[Rcpp::export]]
List long_derivs(const vec& b, const mat& X, const vec& Y, const mat& Z, const vec& beta,
                 const double& shape, const mat& design){
  vec score = get_long_score(b, X, Y, Z, beta, shape, design);
  mat Hessian = get_long_hess(b, X, Y, Z, beta, shape, design);
  return List::create(_["grad"] = score, _["Hessian"] = Hessian);
}

// First derivative of joint density wrt b.
// [[Rcpp::export]]
vec joint_density_ddb(vec& b, mat& X, vec& Y, mat& Z, vec& beta, double shape, 
                      mat& D, rowvec& S, mat& SS, rowvec& Fi, mat& Fu, double l0i, rowvec& haz,
                      int Delta, double gamma, vec& zeta){
  vec grad = get_long_score(b, X, Y, Z, beta, shape, Z);
  return -grad + -1.0 * (-D.i() * b + Delta * (Fi.t() * gamma) - 
                         gamma * (Fu.t() * (haz.t() % exp(SS * zeta + Fu * (gamma * b)))));
}

// Parameter Updates ------------------------------------------------------
// Shape
// [[Rcpp::export]]
List shape_update(const vec& b, const mat& X, const vec& Y, const mat& Z,
                  const vec& beta, const double& shape, const vec& w, const vec& v,
                  const vec& tau){
  double dig = R::psigamma(shape, 0.), psi1 = R::psigamma(shape, 1.); // 1st & 2nd lgamma derivs
  double H = 1./shape - psi1;  // Hessian not dependent on REs.
  // Score, which we need to take w/ quadrature.
  int mi = Y.size(), gh = w.size();
  vec eta = X * beta + Z * b;
  double score_lhs = sum(log(Y) - dig - eta + 1. + log(shape)), score_rhs = 0.;
  for(int l = 0; l < gh; l++){
    vec eta_l = eta + v[l] * tau;
    score_rhs += w[l] * sum(Y/exp(eta_l));
  }
  double Score = sum(score_lhs - score_rhs);
  return List::create(_["Score"] = Score, _["Hessian"] = H);
}

// [[Rcpp::export]]
List shape_update2(const vec& b, const mat& X, const vec& Y, const mat& Z,
                   const vec& beta, const double& shape, const vec& w, const vec& v,
                   const vec& tau){
  double dig = R::psigamma(shape, 0.), psi1 = R::psigamma(shape, 1.); // 1st & 2nd lgamma derivs
  double H = 1./shape - psi1;  // Hessian not dependent on REs.
  // Score, which we need to take w/ quadrature.
  vec eta = X * beta + Z * b;
  int mi = Y.size(), gh = w.size();
  double score = 0.;
  for(int l = 0; l < gh; l++){
    vec eta_l = eta + v[l] * tau;
    vec mu = exp(eta_l);
    score += w[l] * sum(
      log(Y) - dig - eta_l + 1. + log(shape) - Y/mu 
    );
  }
  return List::create(_["Score"] = score, _["Hessian"] = H);
}

// Survival pair (gamma, zeta) ----------------------------

// Define the conditional expectation
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

// Score AND Hessian via forward differencing.
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

// Baseline hazard \lambda_0 ------------------------
// [[Rcpp::export]]
mat lambdaUpdate(List survtimes, mat& ft, double gamma, vec& zeta,
                 List S, List Sigma, List b, vec& w, vec& v){
  int gh = w.size(), id = b.size();
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

/* *****
 * End *
 * *****/


