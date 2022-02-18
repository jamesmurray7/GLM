#include <RcppArmadillo.h>
#include <math.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
double negbin_ll(vec& Y, vec& mu, double theta){
	vec out = vec(mu.size());
	for(int i = 0; i < mu.size(); i++){
		out[i] = R::dnbinom_mu(Y[i], theta, mu[i], 1);
	}
	return sum(out);
}

// dl(eta)/deta for negative binomial
//[[Rcpp::export]]
vec Score_eta(vec& b, mat& X, vec& Y, mat& Z, vec& beta, double theta){
	vec eta = X * beta + Z * b;
	return (theta * (Y - exp(eta)))/(theta + exp(eta));
}

// Negative log-likelihood for joint density f(Y,T|b,...)
// [[Rcpp::export]]
double joint2_density(vec& b, mat& X, vec& Y, mat& Z, vec& beta, double theta, mat& D,
                     int Delta, rowvec& K, rowvec Fi, double l0i, 
                     mat& KK, mat& Fu, rowvec& haz, double gamma, vec& eta){
	vec this_eta = X * beta + Z * b;
	vec mu = exp(this_eta);
	mu.replace(datum::inf, 1e100); // failsafe
	double temp = 0.0;
	if(Delta == 1) temp = log(l0i);
	// The log likelihood parts
	double nb = negbin_ll(Y, mu, theta);
	double RE = -b.size()/2.0 * log(2.0 * M_PI) - 0.5 * log(det(D)) - 0.5 * as_scalar(b.t() * D.i() * b);
	double surv = temp + as_scalar(Delta * (K * eta + Fi * (gamma * b)) - haz * exp(KK * eta + Fu * (gamma * b)));
	return -1.0 * (sum(nb) + RE + surv);
}

// First derivative of negative log likelihood of joint density w.r.t b
// [[Rcpp::export]]
vec joint_density_ddb(vec& b, mat& X, vec& Y, mat& Z, vec& beta, double theta, mat& D,
                     int Delta, rowvec& K, rowvec Fi, double l0i, 
                     mat& KK, mat& Fu, rowvec& haz, double gamma, vec& eta){
	vec Score_nb_b = Z.t() * Score_eta(b, X, Y ,Z, beta, theta);
	return -1.0* (Score_nb_b -D.i() * b + Delta * gamma * Fi.t() - 
	                      gamma * Fu.t() * (haz.t() % exp(KK * eta + gamma * Fu * b)));
}

// Hessian of negative log likelihood of joint density w.r.t b
// done via forward differencing
// [[Rcpp::export]]
mat joint_density_sdb(vec& b, mat& X, vec& Y, mat& Z, vec& beta, double theta, mat& D,
                     int Delta, rowvec& K, rowvec Fi, double l0i, 
                     mat& KK, mat& Fu, rowvec& haz, double gamma, vec& eta, double eps){
	int n = b.size();
	mat out = zeros<mat>(n, n);
	vec f0 = joint_density_ddb(b, X, Y, Z, beta, theta, D, Delta, K, Fi, l0i, KK, Fu,
	                      haz, gamma, eta);
    for(int i = 0; i < n; i++){
		vec bb = b;
		double xi = std::max(b[i], 1.0);
		bb[i] = b[i] + (eps * xi);
		vec fdiff = joint_density_ddb(bb, X, Y, Z, beta, theta, D, Delta, K, Fi, l0i,
	                           KK, Fu, haz, gamma, eta) - f0;
		out.col(i) = fdiff/(bb[i]-b[i]);
	}
	return 0.5 * (out + out.t());
}

// Function for the Score of \beta
// [[Rcpp::export]]
vec Sbeta(vec& beta, mat& X, vec& Y, mat& Z, vec& b, double theta){
	return X.t() * Score_eta(b, X, Y, Z, beta, theta);
}

// Function for Hessian of \beta
// [[Rcpp::export]]
mat Hbeta(vec& beta, mat& X, vec& Y, mat& Z, vec& b, double theta, double eps){
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




