// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <math.h>
#include <float.h>
using namespace Rcpp;
using namespace arma;

// As it appears in glmmTMB
// double ll_genpois(vec& mu, vec& phi, vec& Y){
// 	vec lam = mu / sqrt(phi);
// 	vec theta = 1. - 1./sqrt(phi);
// 	vec out = log(lam) + (Y - 1.) % log(lam + theta % Y) - lam - theta % Y - lgamma(Y + 1.);
// 	return sum(out);
// }


// Log-likelihoods.
// GP1 from Zamani & Ismail (2012)
// https://www.tandfonline.com/doi/pdf/10.1080/03610926.2011.564742?needAccess=true
// [[Rcpp::export]]
double ll_genpois(vec& mu, double & phi, vec& Y){
	vec frac = (mu + phi * Y) / (1. + phi);
  vec out = log(mu) + (Y - 1.) % log(mu + phi * Y) - lgamma(Y + 1.) - Y * log(1. + phi) - frac;
  return sum(out);
}
 
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
double joint_density(vec& b, mat& X, vec& Y, mat& Z,vec& beta, double & phi, 
                     mat& D, rowvec& S, mat& SS, rowvec& Fi, mat& Fu, double l0i, rowvec& haz,
                     int Delta, double gamma, vec& zeta){
  // Define mu, nu
  vec mu = exp(X * beta + Z * b);
  
  // Consituent parts.
  double ll_GP =   ll_genpois(mu, phi, Y);
  double ll_surv = logfti(b, S, SS, Fi, Fu, l0i, haz, Delta, gamma, zeta);
  
  int q = b.size();
  return -1. * (ll_GP + as_scalar(-(double)q/2.0 * log2pi - 0.5 * log(det(D)) - 0.5 * b.t() * D.i() * b) + ll_surv);
}

// Derivatives ------------------------------------------------------------
vec get_long_score(vec& b, mat& X, vec& Y, mat& Z,vec& beta, double & phi, mat& design){
  // Exploiting that d/dbeta===d/db besides the design measure.
  int q = design.n_cols;
  vec mu = exp(X * beta + Z * b), grad = vec(q);
  for(int qq = 0; qq < q; qq++){
    vec x = design.col(qq);
    grad[qq] = sum(x + (Y - 1.) % (x % mu)/(mu + phi * Y) - x % mu / (phi + 1.));
  }
  return grad;
}

mat get_long_hess(vec& b, mat& X, vec& Y, mat& Z,vec& beta, double & phi, mat& design){
  int mi = design.n_rows, q = design.n_cols;
  vec mu = exp(X * beta + Z * b);
  mat H = zeros<mat>(q, q);
  vec inner_part = square(mu + phi * Y);
  for(int j = 0; j < mi; j++){
    rowvec xjT = design.row(j);
    vec xj = xjT.t();
    H += (phi * (Y[j] - 1.) * Y[j] / (inner_part[j]) - mu[j] / (phi + 1.)) * xj * xjT;
  }
  return H;
}

// [[Rcpp::export]]
List long_derivs(vec& b, mat& X, vec& Y, mat& Z,vec& beta, double & phi, mat& design){
  vec score = get_long_score(b,X,Y,Z,beta,phi,design);
  mat Hessian = get_long_hess(b,X,Y,Z,beta,phi,design);
  return List::create(_["grad"] = score, _["Hessian"] = Hessian);
}

// First and second derivatives of joint density wrt b.
// [[Rcpp::export]]
vec joint_density_ddb(vec& b, mat& X, vec& Y, mat& Z, vec& beta, double & phi, 
                     mat& D, rowvec& S, mat& SS, rowvec& Fi, mat& Fu, double l0i, rowvec& haz,
                     int Delta, double gamma, vec& zeta){
  // int q = Z.n_cols;
  vec mu = exp(X * beta + Z * b);//, grad = vec(q);
  // for(int qq = 0; qq < q; qq++){
  //   vec Zq = Z.col(qq);
  //   grad[qq] = sum(Zq + (Y - 1.) % (Zq % mu)/(mu + phi * Y) - Zq % mu / (phi + 1.));
  // }
  vec grad = get_long_score(b, X, Y, Z, beta, phi, Z);
  
  return -grad + -1.0 * (-D.i() * b + Delta * (Fi.t() * gamma) - 
                         gamma * (Fu.t() * (haz.t() % exp(SS * zeta + Fu * (gamma * b)))));
}

// [[Rcpp::export]]
vec joint_density_sdb(vec& b, mat& X, vec& Y, mat& Z,vec& beta, double & phi, 
                      mat& D, rowvec& S, mat& SS, rowvec& Fi, mat& Fu, double l0i, rowvec& haz,
                      int Delta, double gamma, vec& zeta){
  // int mi = Z.n_rows, q = Z.n_cols;
  // vec mu = exp(X * beta + Z * b);
  // mat H = zeros<mat>(q, q);
  // vec inner_part = square(mu + phi * Y);
  // for(int j = 0; j < mi; j++){
  //   rowvec ZjT = Z.row(j);
  //   vec Zj = ZjT.t();
  //   H += (phi * (Y[j] - 1.) * Y[j] / (inner_part[j]) - mu[j] / (phi + 1.)) * Zj * ZjT;
  // }
  mat H = get_long_hess(b, X, Y, Z, beta, phi, Z);
  return -H + -1. * (-D.i() - pow(gamma, 2.) * (diagmat(haz.t() % exp(SS * zeta + Fu * (gamma * b))) * Fu).t() * Fu);
} 

// Parameter updates ------------------------------------------------------
// Dispersion, phi
// [[Rcpp::export]]
List phi_update(vec& b, mat& X, vec& Y, mat& Z,vec& beta, double & phi,
                vec& w, vec& v, vec& tau){
  int gh = w.size();
  double rhs = sum(Y)/(1.+phi), Score = 0., Hess = 0.;
  vec eta = X * beta + Z * b;
  for(int l = 0; l < gh; l++){
    vec eta_l = eta + tau * v[l];
    vec mu = exp(eta_l);
    Score += w[l] * sum(
      ((Y - 1.) % Y)/(mu + phi * Y) + (mu - Y)/(pow(phi + 1., 2.))
    );
    Hess += w[l] * sum(
      (2. * (Y - mu))/(pow(phi + 1., 3.)) - (square(Y) % (Y - 1.))/(square(mu + phi * Y))
    );
  }
  return List::create(_["Score"] = Score - rhs, _["Hessian"] = rhs + Hess);
}
