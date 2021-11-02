#include <RcppArmadillo.h>
#include <math.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// Joint Likelihood ----------------------------------------------------
// [[Rcpp::export]]
double ll(vec& b, const colvec& Y, const colvec& lfactY, const mat& X, const mat& Z,
    		  const mat& D, const vec& theta,   
    		  const rowvec& K, const int Delta, const double l0i, const rowvec& Fi,
    		  const rowvec& l0u, const mat& Fu, const vec& beta, const vec& eta,
    		  const vec& gr, const rowvec& rvFi, const int nK, const int q){
	double temp = 0.0;
	if(Delta == 1) temp = log(l0i);
	return -1.0 * as_scalar(
	   sum(lgamma(theta + Y) - lgamma(theta) - lfactY) + Y.t() * (X * beta + Z * b) + 
		 theta.t() * log(theta) - (theta+Y).t() * (log(exp(X*beta + Z*b)+theta))+ 
	    - q * log(2.0 * M_PI) - 0.5 * log(det(D)) - 0.5 * b.t() * D.i() * b + 
	    temp + Delta * (K * eta + rvFi * (gr % b)) - l0u * (kron(exp(K * eta), exp(repmat(Fu, 1, nK) * (gr % b))))
	 );
}

// Gradient function - NB some redundant arguments are included to allow for ease of use in ucminf
// [[Rcpp::export]]
colvec gradll(vec& b, const colvec& Y, const colvec& lfactY, const mat& X, const mat& Z,
              const mat& D, const vec& theta,   
              const rowvec& K, const int Delta, const double l0i, const rowvec& Fi,
              const rowvec& l0u, const mat& Fu, const vec& beta, const vec& eta,
              const vec& gr, const rowvec& rvFi, const int nK, const int q){
	return -1.0 * (Z.t() * Y - Z.t() * ((Y + theta) % exp(X * beta + Z * b) / (exp(X * beta + Z * b) + theta)) 
	               - D.i() * b + Delta * rvFi.t() % gr + 
	               - repmat(Fu, 1, nK).t() * (l0u.t() % (kron(exp(K * eta), exp(repmat(Fu, 1, nK) * (gr % b))))) % gr);

}

// Just functions for the negative binomial sub-model ----------------------------
// Likelihood
// [[Rcpp::export]]
double nb_ll(vec& b, const colvec& Y, const colvec& lfactY, const mat& X, const mat& Z,
			 const mat& D, const vec& beta, const vec& theta, const int q){
	
	return -1.0 * as_scalar(
		sum(lgamma(theta + Y) - lgamma(theta) - lfactY) + Y.t() * (X * beta + Z * b) + 
		theta.t() * log(theta) - (theta+Y).t() * (log(exp(X*beta + Z*b)+theta)) + 
		-q/2.0 * log(2.0 * M_PI) - 0.5 * log(det(D)) - 0.5 * b.t() * D.i() * b
	);
}

// negative binomial first derivative wrt bi, note redundant arguments for sake of ease to call with ucminf
// [[Rcpp::export]]
colvec nb_grad(vec& b, const colvec& Y, const colvec& lfactY, const mat& X, const mat& Z,
			 const mat& D, const vec& beta, const vec& theta, const int q){
	return -1.0 * (Z.t() * Y - Z.t() * ((Y + theta) % exp(X * beta + Z * b) / (exp(X * beta + Z * b) + theta)) - D.i() * b);
}

// negative binomial second derivate wrt bi
// [[Rcpp::export]]
mat nb_hess(vec& b, const mat& X, const mat& Z, const colvec& Y, const mat& D, const vec& beta,
				const vec& theta){
    vec mu = exp(X * beta + Z * b);
    vec t1 = (Y + theta) % mu;
    vec t2 = mu + theta;
    mat dmat1 = diagmat(t1 / t2);
    mat dmat2 = diagmat((t1 % mu)/(t2%t2));
    
	return -1.0 * ( (dmat1 * Z).t() * Z - (dmat2 * Z).t() * Z ) - D.i() ;
}
  
// ------
// Joint likelihood: Hessian
// [[Rcpp::export]]
mat hessll(vec& b, const mat& X, const mat& Z, const colvec& Y, const mat& D,
		       const vec& theta, const rowvec& K, const vec& l0u, const mat& Fu, 
           const vec& beta, const vec& eta, const vec& gr, const int nK){
  mat nb = nb_hess(b, X, Z, Y, D, beta, theta);
  mat diaggr = diagmat(gr);
  vec kernel = l0u % kron(exp(K * eta), exp(repmat(Fu, 1, nK) * (gr % b)));
  mat diagkern = diagmat(kernel);
  return nb + (-diaggr * repmat(Fu, 1, nK).t() * diagkern * repmat(Fu, 1, nK) * diaggr);
}


// Version with split-out _k's for theta (to be found via numDeriv grad and hess) ( not currently working :|)
// [[Rcpp::export]]
double nb_ll_scalar_theta(double theta, const colvec& Y, const colvec& lfactY, const mat& X, const mat& Z,
			 const mat& D, const vec& beta, const vec& b, const int q){
	
	return -1.0 * as_scalar(
		sum(lgamma(theta + Y) - lgamma(theta) - lfactY) + Y.t() * (X * beta + Z * b) + 
		theta * log(theta) - (theta+Y).t() * (log(exp(X*beta + Z*b)+theta)) + 
		-q/2.0 * log(2.0 * M_PI) - 0.5 * log(det(D)) - 0.5 * b.t() * D.i() * b
	);
}












