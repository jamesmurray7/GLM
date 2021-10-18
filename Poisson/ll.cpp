#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double ll(vec& b, const colvec& Y, const colvec& lfactY, const mat& X, const mat& Z, const mat& D,
		  const rowvec& K, const int Delta, const double l0i, const rowvec& Fi,
		  const rowvec& l0u, const mat& Fu, const vec& beta, const vec& eta,
		  const vec& gr, const rowvec& rvFi, const int nK, const int q){
	double temp = 0.0;
	if(Delta == 1) temp = log(l0i);
	return -1.0 * as_scalar(sum(-lfactY) + sum(-exp(X * beta + Z * b)) + Y.t() * (X * beta + Z * b) + 
	                 - q * log(2.0 * M_PI) - 0.5 * log(det(D)) - 0.5 * b.t() * D.i() * b + 
	                 temp + Delta * (K * eta + rvFi * (gr % b)) - l0u * (kron(exp(K * eta), exp(repmat(Fu, 1, nK) * (gr % b)))));
}

// Gradient function
// [[Rcpp::export]]
colvec gradll(vec& b, const colvec& Y, const mat& X, const mat& Z, const mat& D,
		     const rowvec& K, const int Delta, const double l0i,
		     const rowvec& l0u, const mat& Fu, const vec& beta, const vec& eta,
		     const vec& gr, const rowvec& rvFi, const int nK){
	return -1.0 * (-Z.t() * exp(X * beta + Z * b) + Z.t() * Y - D.i() * b + Delta * rvFi.t() % gr + 
	                 - repmat(Fu, 1, nK).t() * (l0u.t() % (kron(exp(K * eta), exp(repmat(Fu, 1, nK) * (gr % b))))) % gr);

}
