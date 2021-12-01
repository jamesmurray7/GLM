#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


// 1. Function for log density via R::dbinom. 
// This is only marginally faster than R's implementation, but included for ease of use with later functions
// [[Rcpp::export]]
double log_density(vec& eta, vec& Y){
	vec mu = exp(eta)/(1+exp(eta));
	vec out = vec(mu.size());
	for(int i = 0; i < eta.size(); i++){
		out[i] = R::dbinom(Y[i],1,mu[i],1);
	}
	return sum(out);
}

// 2. Function for the Score of \beta
// [[Rcpp::export]]
vec Sbeta(vec& beta, mat& X, vec& Y, mat& Z, vec& b){
	vec eta = X * beta + Z * b;
	return X.t() * (Y - exp(eta)/(exp(eta) + 1));
}

// 3. Function for the Hessian of \beta
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

// 4. Function for ll of full joint density of f(Y,T|X,b,...; params).
// [[Rcpp::export]]
double joint_density(vec& b, mat& X, vec& Y, mat& Z, vec& beta, mat& D,           // Longit binary + REs
					 int Delta, rowvec& K, rowvec& Fi, double l0i,                // Survival
					 mat& KK, mat& Fu, rowvec& haz, double gamma, vec& eta){      // -"-
	vec this_eta = X * beta + Z * b;
	double bin = log_density(this_eta, Y);
	int q = b.size();
	double temp = 0.0;
	if(Delta == 1) temp = log(l0i);
	return -bin + -1.0*as_scalar(-q/2.0 * log(2.0 * M_PI) - 0.5 * log(det(D)) - 0.5 * b.t() * D.i() * b + 
	                       temp + Delta * (K * eta + Fi * (gamma * b)) - haz * exp(KK * eta + Fu * (gamma * b)));
}

// 5. Function for first derivative of logf(Y,T|...) w.r.t b.
// [[Rcpp::export]]
vec joint_density_ddb(vec& b, mat& X, vec& Y, mat& Z, vec& beta, mat& D,           // Longit binary + REs
					 int Delta, rowvec& K, rowvec& Fi, double l0i,                 // Survival
					 mat& KK, mat& Fu, rowvec& haz, double gamma, vec& eta){
	vec Score_bin_b = Z.t() * (Y - exp(X * beta + Z * b)/(exp(X * beta + Z * b) + 1));
	return -Score_bin_b + -1.0 * (-D.i() * b + Delta * gamma * Fi.t() - 
	                      gamma * Fu.t() * (haz.t() % exp(KK * eta + gamma * Fu * b)));
}

// 6. Function for the second derivative of logf(Y,T|...) w.r.t b.
//    This is done via forward differencing.
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


// 7. Function for first derivative of logf(Y,T|...) w.r.t (gamma, eta).
// [[Rcpp::export]]
vec Sgammaeta(vec gammaeta, vec& b, int Delta, rowvec& Fi, 
		      vec& K, mat& KK, mat& Fu, vec& haz, mat& S, 
		      vec& w, vec& v){
	vec out = vec(gammaeta.size());
	int gh = w.size();
	vec eta = gammaeta.tail(2);
	double gamma = as_scalar(gammaeta.head(1));
	vec tau = diagvec(Fu * S * Fu.t());
	uvec zzz = find(tau == 0.0);
	vec tau2 = pow(pow(gamma, 2.0) * tau, -0.5);
	tau2.elem(zzz).zeros();
	double gamma_rhs = 0.0;
	vec eta_rhs = vec(eta.size());
	for(int l = 0; l < gh; l++){
		vec xi = haz % exp(KK * eta + gamma * Fu * b + pow(pow(gamma, 2.0) * tau, 0.5) * v[l]);
		gamma_rhs += as_scalar(w[l] * xi.t() * (Fu * b) + gamma * w[l] * v[l] * (xi % tau2).t() * tau);
		eta_rhs += w[l] * KK.t() * xi;
	}
	out.head(1) = as_scalar(Delta * Fi * b) - gamma_rhs;
	out.tail(2) = Delta * K - eta_rhs;
	return out;
}

// 8. Function for second derivative of logf(Y,T|...) w.r.t (gamma, eta).
//    This is done via forward differencing.
// [[Rcpp::export]]
mat Hgammaeta(vec gammaeta, vec& b, int Delta, rowvec& Fi, 
		      vec& K, mat& KK, mat& Fu, vec& haz, mat& S, 
		      vec& w, vec& v, double eps){
	int n = gammaeta.size();
	mat out = zeros<mat>(n, n);
	vec f0 = Sgammaeta(gammaeta, b, Delta, Fi, K, KK, Fu, haz, S, w, v);
	for (int i = 0; i < n; i++){
		vec ge = gammaeta;
		double xi = std::max(gammaeta[i], 1.0);
		ge[i] = gammaeta[i] + (eps * xi);
		vec fdiff = Sgammaeta(ge, b, Delta, Fi, K, KK, Fu, haz, S, w, v) - f0;
		out.col(i) = fdiff/(ge[i]-gammaeta[i]);
	}
	return 0.5 * (out + out.t());
}

// 9. Function for the update to the baseline hazard, \lambda_0(u)
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

