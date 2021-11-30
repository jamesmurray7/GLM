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

// 4. Function for full joint density of Y,T|X,b,...;params
// [[Rcpp::export]]
double joint_density(vec& b, mat& X, vec& Y, mat& Z, vec& beta, mat& D,           // Longit binary + REs
					 int Delta, rowvec& K, rowvec& Fi, double l0i,                // Survival
					 mat& KK, mat& Fu, rowvec& haz, double gamma, vec& eta){      // -"-
	
}
