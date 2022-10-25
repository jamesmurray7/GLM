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

// Specifically for obtaining pmf in simulations --------------------------
// [[Rcpp::export]]
double GP1_pmf_scalar(double mu, double & phi, double & Y){
  double L = mu * pow(mu + phi* Y, Y - 1.) * exp(-(mu+phi*Y)/(1.+phi)) / (pow(1. + phi, Y) * tgamma(Y + 1.));
  vec LL = vec(1, fill::value(L));
  LL.replace(datum::inf, 1e100);
  return as_scalar(LL);
}

// Log-likelihoods.
// 'GP1' from Zamani & Ismail (2012)
// https://www.tandfonline.com/doi/pdf/10.1080/03610926.2011.564742?needAccess=true
// [[Rcpp::export]]
double ll_genpois(vec& mu, vec& phi, vec& Y){
	vec frac = (mu + phi % Y) / (1. + phi);
  vec out = log(mu) + (Y - 1.) % log(mu + phi % Y) - lgamma(Y + 1.) - Y % log(1. + phi) - frac;
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
double joint_density(vec& b, mat& X, vec& Y, mat& Z, mat& W,
                     vec& beta, vec& phi, 
                     mat& D, rowvec& S, mat& SS, rowvec& Fi, mat& Fu, double l0i, rowvec& haz,
                     int Delta, double gamma, vec& zeta){
  // Define mu, nu
  vec mu = exp(X * beta + Z * b);
  vec phivec = W * phi;
  // Consituent parts.
  double ll_GP =   ll_genpois(mu, phivec, Y);
  double ll_surv = logfti(b, S, SS, Fi, Fu, l0i, haz, Delta, gamma, zeta);
  int q = b.size();
  return -1. * (ll_GP + as_scalar(-(double)q/2.0 * log2pi - 0.5 * log(det(D)) - 0.5 * b.t() * D.i() * b) + ll_surv);
}

// Derivatives ------------------------------------------------------------
vec get_long_score(vec& b, mat& X, vec& Y, mat& Z, vec& beta, vec& phivec, mat& design){
  // Exploiting that d/dbeta===d/db besides the design measure.
  int q = design.n_cols;
  vec mu = exp(X * beta + Z * b), grad = vec(q);
  for(int qq = 0; qq < q; qq++){
    vec x = design.col(qq);
    grad[qq] = sum(x + (Y - 1.) % (x % mu)/(mu + phivec % Y) - x % mu / (phivec + 1.));
  }
  return grad;
}

mat get_long_hess(vec& b, mat& X, vec& Y, mat& Z,vec& beta, vec& phivec, mat& design){
  int mi = design.n_rows, q = design.n_cols;
  vec mu = exp(X * beta + Z * b);
  mat H = zeros<mat>(q, q);
  vec inner_part = square(mu + phivec % Y);
  for(int j = 0; j < mi; j++){
    rowvec xjT = design.row(j);
    vec xj = xjT.t();
    H += (phivec[j] * (Y[j] - 1.) * Y[j] / (inner_part[j]) - mu[j] / (phivec[j] + 1.)) * xj * xjT;
  }
  return H;
}

// [[Rcpp::export]]
List long_derivs(vec& b, mat& X, vec& Y, mat& Z,vec& beta, vec& phivec, mat& design){
  vec score = get_long_score(b, X, Y, Z, beta, phivec, design);
  mat Hessian = get_long_hess(b, X, Y, Z, beta, phivec, design);
  return List::create(_["grad"] = score, _["Hessian"] = Hessian);
}

// First and second derivatives of joint density wrt b.
// [[Rcpp::export]]
vec joint_density_ddb(vec& b, mat& X, vec& Y, mat& Z, mat& W, 
                      vec& beta, vec & phi, 
                      mat& D, rowvec& S, mat& SS, rowvec& Fi, mat& Fu, double l0i, rowvec& haz,
                      int Delta, double gamma, vec& zeta){
  // int q = Z.n_cols;
  vec phivec = W * phi;
  vec grad = get_long_score(b, X, Y, Z, beta, phivec, Z);
  
  return -grad + -1.0 * (-D.i() * b + Delta * (Fi.t() * gamma) - 
                         gamma * (Fu.t() * (haz.t() % exp(SS * zeta + Fu * (gamma * b)))));
}

// [[Rcpp::export]]
vec joint_density_sdb(vec& b, mat& X, vec& Y, mat& Z, mat& W,             // Not really used!
                      vec& beta, vec & phi,  
                      mat& D, rowvec& S, mat& SS, rowvec& Fi, mat& Fu, double l0i, rowvec& haz,
                      int Delta, double gamma, vec& zeta){
  vec phivec = W * phi;
  mat H = get_long_hess(b, X, Y, Z, beta, phivec, Z);
  return -H + -1. * (-D.i() - pow(gamma, 2.) * (diagmat(haz.t() % exp(SS * zeta + Fu * (gamma * b))) * Fu).t() * Fu);
} 

// Parameter updates ------------------------------------------------------
// Dispersion, phi
// [[Rcpp::export]]
List phi_update(vec& b, mat& X, vec& Y, mat& Z, mat& W, 
                vec& beta, vec & phi,
                vec& w, vec& v, vec& tau){
  int gh = w.size(), ww = W.n_cols, mi = W.n_rows;
  // Rcout << "Matrix W is " << mi << " x " << ww << std::endl;
  vec eta = X * beta + Z * b;
  if(ww == 1){ // Special case, which we can do a little more quickly...
    // Rcout << "Doing one-dimensional" << std::endl;
    double phi_ = phi.at(0);
    double rhs = sum(Y)/(1.+phi_), lhs = sum(Y)/(pow(1.+phi_,2.)), Score = 0., Hess = 0.;
    for(int l = 0; l < gh; l++){
      vec eta_l = eta + tau * v[l];
      vec mu = exp(eta_l);
      Score += w[l] * sum(
        ((Y - 1.) % Y)/(mu + phi_ * Y) + (mu - Y)/(pow(phi_ + 1., 2.))
      );
      Hess += w[l] * sum(
        (2. * (Y - mu))/(pow(phi_ + 1., 3.)) - (square(Y) % (Y - 1.))/(square(mu + phi_ * Y))
      );
    }
    return List::create(_["Score"] = Score - rhs, _["Hessian"] = lhs + Hess);
  }else{ // Else need to appraise each column and each row.
    // Rcout << "Doing w-dimensional" << std::endl;
    vec phivec = W * phi; //, Score = vec(ww);
    mat Hess = zeros<mat>(ww, ww), Score = mat(mi, ww);
    // Loop over quadrature nodes
    for(int l = 0; l < gh; l++){
      vec eta_l = eta + tau * v[l];
      vec mu = exp(eta_l);
      // Score
      for(int c = 0; c < ww; c++){
        vec Wcol = W.col(c);
        Score.col(c) += w[l] * (                          // Slight inefficiency here!
          Wcol % (Y - 1.) % Y / (mu + phivec % Y) - Wcol % Y / (1. + phivec) + 
            Wcol % (mu - Y) / square(1. + phivec)
        );
      }
      // Rcout << "Score: " << Score << std::endl;
      // Hessian
      for(int j = 0; j < mi; j++){
        rowvec WjT = W.row(j);
        vec Wj = WjT.t();
        Hess += w[l] * (
          (Y[j] / pow(1.+phivec[j], 2.) + 2. * (Y[j] - mu[j]) / pow(1. + phivec[j], 3.) -
            (Y[j] - 1.) * pow(Y[j], 2.) / pow(mu[j] + phivec[j]* Y[j], 2.)) * Wj * WjT
        );
      }
    }
    return List::create(_["Score"] = sum(Score, 0), _["Hessian"] = Hess);
  }
  
  
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
