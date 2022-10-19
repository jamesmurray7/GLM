// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <math.h>
#include <float.h>
using namespace Rcpp;
using namespace arma;

// Specifically for obtaining pmf in simulations --------------------------
/*
 f(y|\theta,\lambda) =
 \frac{\theta(theta+\lambda y)^{y-1}e^{-\theta-\lambda y}}{y \!}
 */
// [[Rcpp::export]]
double GP_pmf_scalar(double mu, double & phi, double & Y){
  double theta = mu/sqrt(phi);
  double lambda = 1. - 1./sqrt(phi);
  double L = theta * pow(theta + lambda * Y, Y - 1.) * exp(-theta - lambda * Y)/tgamma(Y+1);
  vec LL = vec(1, fill::value(L));
  LL.replace(datum::inf, 1e100);
  return as_scalar(LL);
}

// Log-likelihoods.
// GP from glmmTMB (https://github.com/glmmTMB/glmmTMB/blob/master/glmmTMB/src/distrib.h).
// I also replace {theta -> lambda, lambda -> nu} hereafter.
// [[Rcpp::export]]
double ll_genpois(vec& mu, double & phi, vec& Y){
  int mi = Y.size();
  vec lambda = mu/sqrt(phi);
  vec nu = vec(mi, fill::value(1. - 1./sqrt(phi)));
	vec out = log(lambda) + (Y - 1.) % log(lambda + nu % Y) - lambda - nu % Y - lgamma(Y + 1.);
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
vec get_long_score(vec& b, mat& X, vec& Y, mat& Z, vec& beta, double & phi, mat& design){
  // Exploiting that d/dbeta===d/db besides the design measure.
  // phi = exp(delta).
  int q = design.n_cols;
  double sqrphi = sqrt(phi);
  vec mu = exp(X * beta + Z * b), grad = vec(q), denom = mu - Y + sqrphi * Y;
  for(int qq = 0; qq < q; qq++){
    vec x = design.col(qq);
    grad[qq] = sum(
      x + x % (Y - 1.) % mu / denom - x % mu / sqrphi
    );
  }
  return grad;
}

mat get_long_hess(vec& b, mat& X, vec& Y, mat& Z, vec& beta, double & phi, mat& design){
  int mi = design.n_rows, q = design.n_cols;
  double sqrphi = sqrt(phi);
  vec mu = exp(X * beta + Z * b);
  vec denom = mu - Y + sqrphi * Y, numer = (sqrphi - 1.) * (Y - 1.) % Y % mu;
  mat H = zeros<mat>(q, q);
  for(int j = 0; j < mi; j++){
    rowvec xjT = design.row(j);
    vec xj = xjT.t();
    H += (numer[j]/pow(denom[j], 2.) - mu[j]/sqrphi) * xj * xjT;
  }
  return H;
}

// [[Rcpp::export]]
List long_derivs(vec& b, mat& X, vec& Y, mat& Z,vec& beta, double & phi, mat& design){
  vec score = get_long_score(b,X,Y,Z,beta,phi,design);
  mat Hessian = get_long_hess(b,X,Y,Z,beta,phi,design);
  return List::create(_["grad"] = score, _["Hessian"] = Hessian);
}

// gradient of f(Y|b) wrt b
// [[Rcpp::export]]
vec joint_density_ddb(vec& b, mat& X, vec& Y, mat& Z, vec& beta, double & phi, 
                      mat& D, rowvec& S, mat& SS, rowvec& Fi, mat& Fu, double l0i, rowvec& haz,
                      int Delta, double gamma, vec& zeta){
  // int q = Z.n_cols;
  vec mu = exp(X * beta + Z * b);//, grad = vec(q);
  vec grad = get_long_score(b, X, Y, Z, beta, phi, Z);
  
  return -grad + -1.0 * (-D.i() * b + Delta * (Fi.t() * gamma) - 
                         gamma * (Fu.t() * (haz.t() % exp(SS * zeta + Fu * (gamma * b)))));
}

// Think this is (too) accurate...
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
// Dispersion, delta
// [[Rcpp::export]]
List delta_update(vec& b, mat& X, vec& Y, mat& Z, vec& beta, double & delta,
                  vec& w, vec& v, vec& tau){ 
  int gh = w.size();
  double phi = exp(delta), Score = 0., Hess = 0.;
  double sqrphi = sqrt(phi);
  vec eta = X * beta + Z * b;
  for(int l = 0; l < gh; l++){
    vec eta_l = eta + tau * v[l];
    vec mu = exp(eta_l);
    vec numer = (Y - 1.) % (Y - mu);
    vec denom = (sqrphi - 1.) * Y + mu;
    Score += w[l] * sum(
      -0.5 + numer / (2. * denom) + mu/(2. * sqrphi) - Y/(2. * sqrphi)
    );
    Hess += w[l] * sum(
      -sqrphi * numer % Y / (4. * pow(denom, 2.)) - mu/(4. * sqrphi)
    );
  }
  return List::create(_["Score"] = Score,
                      _["Hessian"] = Hess + sum(Y)/(4. * sqrphi));
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
