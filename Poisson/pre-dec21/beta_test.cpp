#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// crossprod(Y, X %*% beta + Z %*% b) - exp(sum(X %*% beta + Z %*% b))

// [[Rcpp::export]]
double beta_ll(vec& beta, mat& X, mat& Z, colvec& Y, vec& b){
	return as_scalar(Y.t() * (X * beta + Z * b) - sum(exp(X * beta + Z * b)));
} 

// [[Rcpp::export]]
double beta_ll_quadrature(vec& beta, mat& X, mat& Z, colvec& Y, vec& b,
							vec& tau, vec& w, vec& v, int gh){
    double rhs = 0.0;
    for(int l = 0; l < gh; l++){
		rhs += as_scalar(w[l] * sum(exp(X * beta + Z * b + tau * v[l])));
	}
	return as_scalar(Y.t() * (X * beta + Z * b) - rhs);
} 
