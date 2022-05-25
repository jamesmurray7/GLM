#include <RcppArmadillo.h>
#include <math.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
double negbin_ll(vec& Y, vec& mu, vec& theta){ // Ensure theta is exponenitated!
	vec out = vec(mu.size());
	for(int i = 0; i < mu.size(); i++){
		out[i] = R::dnbinom_mu(Y[i], theta[i], mu[i], 1);
	}
	return sum(out);
}

// dl(eta)/deta for negative binomial
//[[Rcpp::export]]
vec Score_eta(vec& b, mat& X, vec& Y, mat& Z, vec& beta, vec& theta){ // Ensure theta is exponenitated!
	vec eta = X * beta + Z * b;
	return (theta % (Y - exp(eta)))/(theta + exp(eta));
}

// Negative log-likelihood for joint density f(Y,T|b,...)
// [[Rcpp::export]]
double joint_density(vec& b, mat& X, vec& Y, mat& Z, mat& W, vec& beta, vec& alpha, 
                     mat& D, int Delta, rowvec& S, rowvec Fi, double l0i, 
                     mat& SS, mat& Fu, rowvec& haz, double gamma, vec& zeta){
	vec eta = X * beta + Z * b;
  vec theta = exp(W * alpha);
	vec mu = exp(eta);
	// mu.replace(datum::inf, 1e100); // failsafe
	double temp = 0.0;
	if(Delta == 1) temp = log(l0i);
	// The log likelihood parts
	double nb = negbin_ll(Y, mu, theta);
	double RE = -b.size()/2.0 * log(2.0 * M_PI) - 0.5 * log(det(D)) - 0.5 * as_scalar(b.t() * D.i() * b);
	double surv = temp + as_scalar(Delta * (S * zeta + Fi * (gamma * b)) - haz * exp(SS * zeta + Fu * (gamma * b)));
	return -1.0 * (sum(nb) + RE + surv);
}

// First derivative of negative log likelihood of joint density w.r.t b
// [[Rcpp::export]]
vec joint_density_ddb(vec& b, mat& X, vec& Y, mat& Z, mat& W, vec& beta, vec& alpha, 
                      mat& D, int Delta, rowvec& S, rowvec Fi, double l0i, 
                      mat& SS, mat& Fu, rowvec& haz, double gamma, vec& zeta){
  vec theta = exp(W * alpha);
	vec Score_nb_b = Z.t() * Score_eta(b, X, Y ,Z, beta, theta);
	return -1.0* (Score_nb_b -D.i() * b + Delta * gamma * Fi.t() - 
	                      gamma * Fu.t() * (haz.t() % exp(SS * zeta + gamma * Fu * b)));
}

// Hessian of negative log likelihood of joint density w.r.t b
// done via forward differencing
// [[Rcpp::export]]
mat joint_density_sdb(vec& b, mat& X, vec& Y, mat& Z, mat& W, vec& beta, vec& alpha, 
                      mat& D, int Delta, rowvec& S, rowvec Fi, double l0i, 
                      mat& SS, mat& Fu, rowvec& haz, double gamma, vec& zeta, double eps){
	int n = b.size();
	mat out = zeros<mat>(n, n);
	vec f0 = joint_density_ddb(b, X, Y, Z, W, beta, alpha, D, Delta, S, Fi, l0i, SS, Fu,
	                      haz, gamma, zeta);
    for(int i = 0; i < n; i++){
		vec bb = b;
		double xi = std::max(b[i], 1.0);
		bb[i] = b[i] + (eps * xi);
		vec fdiff = joint_density_ddb(bb, X, Y, Z, W, beta, alpha, D, Delta, S, Fi, l0i, SS, Fu,
                                haz, gamma, zeta) - f0;
		out.col(i) = fdiff/(bb[i]-b[i]);
	}
	return 0.5 * (out + out.t()); // ensure symmtry...
}

// Function for the Score of \beta
// [[Rcpp::export]]
vec Sbeta(vec& beta, List X, List Y, List Z, List b, List theta){
  int n = X.size();
  int p = beta.size();
  vec Score = vec(p);
  for(int i = 0; i < n; i++){
    mat Xi = X[i];
    vec Yi = Y[i];
    mat Zi = Z[i];
    vec bi = b[i];
    vec th = theta[i];
    Score += Xi.t() * Score_eta(bi, Xi, Yi, Zi, beta, th);
  }
	return Score;
}

// Function for Hessian of \beta
// [[Rcpp::export]]
mat Hbeta(vec& beta, List X, List Y, List Z, List b, List theta, double eps){
	int n = beta.size();
	mat out = zeros<mat>(n, n);
	vec f0 = Sbeta(beta, X, Y, Z, b, theta);
	for(int i = 0; i < n; i++){
		vec bb = beta;
		double xi = std::max(beta[i], 1.0);
		bb[i] = beta[i] + (eps * xi);
		vec fdiff = Sbeta(bb, X, Y, Z, b, theta) - f0;
		out.col(i) = fdiff/(bb[i]-beta[i]);
	}
	return 0.5 * (out + out.t());
}

// Update to dispersion parameter, theta
double ll_theta_quad(vec& alpha, vec& beta, mat& X, vec& Y, mat& Z, mat& W, vec& b, mat& S,
                     vec& w, vec& v){
  vec rhs = vec(Y.size());
  vec tau = sqrt(diagvec(Z * S * Z.t()));
  vec eta = X * beta + Z * b;
  vec theta = exp(W * alpha);
  for(int l = 0; l < w.size(); l++){
    vec this_eta = eta + v[l] * tau;
    rhs += w[l] * log(exp(this_eta) + theta);
  }
  vec out = vec(Y.size());
  for(int i = 0; i < Y.size(); i++){
    out[i] = lgamma(theta[i] + Y[i]) - lgamma(theta[i]) + theta[i] * log(theta[i]) - (theta[i] + Y[i]) * rhs[i];
  }
  return sum(out);
}

// [[Rcpp::export]]
vec Salpha(vec& alpha, vec& beta, mat& X, vec& Y, mat& Z, mat& W, vec& b, mat& S,
           vec& w, vec& v, double eps){
  int a = alpha.size();
  vec out = vec(a);
  double f0 = ll_theta_quad(alpha, beta, X, Y, Z, W, b, S, w, v);
  for(int i = 0; i < a; i++){
    vec aa = alpha;
    double xi = std::max(alpha[i], 1.0);
    aa[i] = alpha[i] + eps * xi;
    double feps = ll_theta_quad(aa, beta, X, Y, Z, W, b, S, w, v) - f0;
    out[i] = feps/(aa[i]-alpha[i]);
  }
  return out;
}

// [[Rcpp::export]]
mat Halpha(vec& alpha, vec& beta, mat& X, vec& Y, mat& Z, mat& W, vec& b, mat& S,
           vec& w, vec& v, double eps1, double eps2){
  int a = alpha.size();
  mat out = zeros<mat>(a,a);
  vec f0 = Salpha(alpha, beta, X, Y, Z, W, b, S, w, v, eps1);
  for(int i = 0; i < a; i++){
    vec aa = alpha;
    double xi = std::max(alpha[i], 1.0);
    aa[i] = alpha[i] + eps2 * xi;
    vec feps = Salpha(aa, beta, X, Y, Z, W, b, S, w, v, eps1) - f0;
    out.col(i) = feps/(aa[i]-alpha[i]);
  }
  return 0.5 * (out + out.t());
}

// Functions for derivatives of logf(Y,T|...) w.r.t (gamma, eta).
// Define the conditional expectation and then take Score AND Hessian via forward differencing
// [[Rcpp::export]]
double Egammaeta(vec& gammaeta, vec& b, mat& S,
                 rowvec& K, mat& KK, mat& Fu, rowvec& Fi, vec& haz, int Delta, vec& w, vec& v){
  double g = as_scalar(gammaeta.head(1));
  vec e = gammaeta.tail(2);
  vec tau = pow(g, 2.0) * diagvec(Fu * S * Fu.t());
  double rhs = 0.0;
  for(int l = 0; l < w.size(); l++){
    rhs += w[l] * as_scalar(haz.t() * exp(KK * e + Fu * (g * b) + v[l] * pow(tau, 0.5)));
  }
  return as_scalar(Delta * (K * e + Fi * (g * b)) - rhs);
}
// [[Rcpp::export]]
vec Sgammaeta(vec& gammaeta, vec& b, mat& S,
              rowvec& K, mat& KK, mat& Fu, rowvec& Fi, vec& haz, int Delta, vec& w, vec& v, double eps){
  vec out = vec(gammaeta.size());
  double f0 = Egammaeta(gammaeta, b, S, K, KK, Fu, Fi, haz, Delta, w, v);
  for(int i = 0; i < gammaeta.size(); i++){
    vec ge = gammaeta;
    double xi = std::max(ge[i], 1.0);
    ge[i] = gammaeta[i] + xi * eps;
    double fdiff = Egammaeta(ge, b, S, K, KK, Fu, Fi, haz, Delta, w, v)-f0;
    out[i] = fdiff/(ge[i]-gammaeta[i]);
  }
  return out;
}
// [[Rcpp::export]]
mat Hgammaeta(vec& gammaeta, vec& b, mat& S,
              rowvec& K, mat& KK, mat& Fu, rowvec& Fi, vec& haz, int Delta, vec& w, vec& v, double eps){
  mat out = zeros<mat>(gammaeta.size(), gammaeta.size());
  vec f0 = Sgammaeta(gammaeta, b, S, K, KK, Fu, Fi, haz, Delta, w, v, eps);
  for(int i = 0; i < gammaeta.size(); i++){
    vec ge = gammaeta;
    double xi = std::max(ge[i], 1.0);
    ge[i] = gammaeta[i] + xi * eps;
    vec fdiff = Sgammaeta(ge, b, S, K, KK, Fu, Fi, haz, Delta, w, v, eps) - f0;
    out.col(i) = fdiff/(ge[i]-gammaeta[i]);
  }
  return 0.5 * (out + out.t());
}


// Function for the update to the baseline hazard, \lambda_0(u)
//     (adapted from my implementation in Bernhardt work)...
// [[Rcpp::export]]
mat lambdaUpdate(List survtimes, vec& ft,
				 double gamma, vec eta, List K, List S, List b, 
				 const vec& w, const vec& v){
  int gh = w.size();
  int id = b.size();
	mat store = zeros<mat>(ft.size(), id); // Initialise the matrix
	// Start loop over i subjects
	for(int i = 0; i < id; i++){
		vec survtimes_i = survtimes[i];    // This id's survived time indices   
		mat Si = S[i];
		vec bi = b[i];
		rowvec Ki = K[i];                  // Start loop over subject i's j survived times     
		for(int j = 0; j < survtimes_i.size(); j++){
			rowvec Fst = NumericVector::create(1.0, ft[j]);
			double tau = as_scalar(pow(gamma, 2.0) * Fst * Si * Fst.t());
			vec rhs = NumericVector::create(0.0, 0.0); // intslope hardcode...
			rhs += gamma * bi;			
			double mu = as_scalar(exp(Ki * eta + Fst * rhs));
			// Rcpp::Rcout << "mu = " << mu << std::endl;								  
			for(int l = 0; l < gh; l++){ // Finally loop over gh nodes
				store(j,i) += as_scalar(w[l] * mu * exp(v[l] * sqrt(tau)));
			}
		}
	}
	return store;
}




