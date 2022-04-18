#include <RcppArmadillo.h>
#include <math.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// Log-likelihood for Poisson via R::dpois(..., log = TRUE).
// [[Rcpp::export]]
double pois_ll(vec& b, mat& X, vec& Y, mat& Z, vec& beta){
	vec eta = X * beta + Z * b;
	vec out = vec(Y.size());
	for(int i = 0; i < Y.size(); i++){
		out[i] = R::dpois(Y[i], exp(eta[i]), 1);
	}
	return sum(out);
}

// conditional expectation wrt REs for Poisson log-likelihood via quadature
vec lfact(vec& v){
	vec out = vec(v.size());
	for(int i = 0; i < v.size(); i++){
		out[i] = lgamma(v[i] + 1.0);
	}
	return out;
}

// [[Rcpp::export]]
double pois_ll_quad(vec& beta, mat& X, vec& Y, mat& Z, vec& b, mat& S,
                    vec& w, vec& v){
	int gh = w.size();
	vec out = vec(Y.size()), rhs = out;
	vec eta = X * beta + Z * b;
	vec tau = sqrt(diagvec(Z * S * Z.t()));
	for(int l = 0; l < gh; l++){
		vec this_eta = eta + v[l] * tau;
		rhs += w[l] * exp(this_eta);
	}
	out = Y % log(rhs) - rhs - lfact(Y);
	return sum(out);


// [[Rcpp::export]]
vec Sbeta_quad(vec& beta, mat& X, vec& Y, mat& Z, vec& b, vec& tau, vec& w, vec& v, long double eps){
	int n = beta.size();
	vec out = vec(n);
	double f0 = pois_ll_quad(beta, X, Y, Z, b, tau, w, v);
	for(int i = 0; i < n; i++){
		vec bb = beta;
		double xi = std::max(1.0, beta[i]);
		bb[i] = beta[i] + (xi * eps);
		double fdiff = pois_ll_quad(bb, X, Y, Z, b, tau, w, v) - f0;
		out[i] = fdiff/(bb[i]-beta[i]);
	}
	return out;
}

// [[Rcpp::export]]
mat Hbeta_quad(vec& beta, mat& X, vec& Y, mat& Z, vec& b, vec& tau, vec& w, vec& v, long double eps){
	int n = beta.size();
	mat out = zeros<mat>(n,n);
	vec f0 = Sbeta_quad(beta, X, Y, Z, b, tau, w, v, eps);
	for(int i = 0; i < n; i++){
		vec bb = beta;
		double xi = std::max(1.0, beta[i]);
		bb[i] = beta[i] + (xi * eps);
		vec fdiff = Sbeta_quad(bb, X, Y, Z, b, tau, w, v, eps) - f0;
		out.col(i) = fdiff/(bb[i]-beta[i]);
	}
	return out;
}

// Score for linear predictor in poisson ll.
// [[Rcpp::export]]
vec Score_eta(vec& b, mat& X, vec& Y, mat& Z, vec& beta){
	vec eta = X * beta + Z * b;
	return Y - exp(eta);
}

// Score of \beta
// [[Rcpp::export]]
vec Sbeta(vec& beta, mat& X, vec& Y, mat& Z, vec& b){
	vec eta = X * beta + Z * b;
	return X.t() * Score_eta(b, X, Y, Z, beta); // watch b<->beta here!
}

// Hessian of \beta, done via forward differencing...
// [[Rcpp::export]]
mat Hbeta(vec& beta, mat& X, vec& Y, mat& Z, vec& b, double eps){
	int n = beta.size();
	mat out = zeros<mat>(n, n);
	vec f0 = Sbeta(beta, X, Y, Z, b);
	for(int i = 0; i < n; i++){
		vec bb = beta;
		double xi = std::max(beta[i], 1.0);
		bb[i] = beta[i] + (eps * xi);
		vec fdiff = Sbeta(bb, X, Y, Z, b) - f0;
		out.col(i) = fdiff/(bb[i]-beta[i]);
	}
	return 0.5 * (out + out.t());
}

// Joint density (negative) log likelihood
// [[Rcpp::export]]
double joint_density(vec& b, mat& X, vec& Y, mat& Z, vec& beta, mat& D,           // Longit binary + REs
					 int Delta, rowvec& K, rowvec& Fi, double l0i,                // Survival
					 mat& KK, mat& Fu, rowvec& haz, double gamma, vec& eta){      // -"-
	double po = pois_ll(b, X, Y, Z, beta);
	int q = b.size();
	double temp = 0.0;
	if(Delta == 1) temp = log(l0i);
	return -po + -1.0 * as_scalar(-q/2.0 * log(2.0 * M_PI) - 0.5 * log(det(D)) - 0.5 * b.t() * D.i() * b + 
	                         temp + Delta * (K * eta + Fi * (gamma * b)) - haz * exp(KK * eta + Fu * (gamma * b)));
}

// First derivative of log likelihood w.r.t b
// [[Rcpp::export]]
vec joint_density_ddb(vec& b, mat& X, vec& Y, mat& Z, vec& beta, mat& D, 
					  int Delta, rowvec& K, rowvec& Fi, double l0i,  
					  mat& KK, mat& Fu, rowvec& haz, double gamma, vec& eta){
	vec Score_po_b = Z.t() * Score_eta(b, X, Y, Z, beta);
	return -Score_po_b + -1.0 * (-D.i() * b + Delta * gamma * Fi.t() - 
	                      gamma * Fu.t() * (haz.t() % exp(KK * eta + gamma * Fu * b)));
}

// Second derivative of log likelihood w.r.t b, done via forward differencing
// [[Rcpp::export]]
mat joint_density_sdb(vec& b, mat& X, vec& Y, mat& Z, vec& beta, mat& D,           // Longit binary + REs
					  int Delta, rowvec& K, rowvec& Fi, double l0i,                 // Survival
					  mat& KK, mat& Fu, rowvec& haz, double gamma, vec& eta, double eps){
	int n = b.size();
	mat out = zeros<mat>(n, n);
	vec f0 = joint_density_ddb(b, X, Y, Z, beta, D, Delta, K, Fi, l0i,
	                           KK, Fu, haz, gamma, eta);
	for(int i = 0; i < n; i++){
		vec bb = b;
		double xi = std::max(b[i], 1.0);
		bb[i] = b[i] + (eps * xi);
		vec fdiff = joint_density_ddb(bb, X, Y, Z, beta, D, Delta, K, Fi, l0i,
	                           KK, Fu, haz, gamma, eta) - f0;
		out.col(i) = fdiff/(bb[i]-b[i]);
	}
	return 0.5 * (out + out.t());
}

// 7. Updates to (gamma, eta)
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



