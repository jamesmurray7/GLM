#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// Joint Likelihood ----------------------------------------------------
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

// Gradient function - NB some redundant arguments are included to allow for ease of use in ucminf
// [[Rcpp::export]]
colvec gradll(vec& b, const colvec& Y, const colvec& lfactY, const mat& X, const mat& Z, const mat& D,
		  const rowvec& K, const int Delta, const double l0i, const rowvec& Fi,
		  const rowvec& l0u, const mat& Fu, const vec& beta, const vec& eta,
		  const vec& gr, const rowvec& rvFi, const int nK, const int q){
	return -1.0 * (-Z.t() * exp(X * beta + Z * b) + Z.t() * Y - D.i() * b + Delta * rvFi.t() % gr + 
	                 - repmat(Fu, 1, nK).t() * (l0u.t() % (kron(exp(K * eta), exp(repmat(Fu, 1, nK) * (gr % b))))) % gr);

}

// Just functions for the poisson sub-model ----------------------------
// Likelihood
// [[Rcpp::export]]
double po_ll(vec& b, const colvec& Y, const colvec& lfactY, const mat& X, const mat& Z,
				const mat& D, const vec& beta, const int q){
	return -1.0 * as_scalar(
		sum(-lfactY) + sum(-exp(X * beta + Z * b)) + Y.t() * (X * beta + Z * b) + 
		-q/2.0 * log(2.0 * M_PI) - 0.5 * log(det(D)) - 0.5 * b.t() * D.i() * b
	);
}

// Poisson first derivative wrt bi
// [[Rcpp::export]]
colvec po_grad(vec& b, const colvec& Y, const colvec& lfactY, const mat& X, const mat& Z,
				const mat& D, const vec& beta, const int q){
	return -1.0 * (-Z.t() * exp(X * beta + Z * b) + Z.t() * Y - D.i() * b);
}

// Poisson second derivate wrt bi
// [[Rcpp::export]]
mat po_hess(vec& b, const mat& X, const mat& Z, const mat& D, const vec& beta){
	vec temp = exp(X * beta + Z * b);
	mat dmat = diagmat(temp);
	return -1.0 * ((dmat * Z).t() * Z) - D.i();
}
  
// ------
// Joint likelihood: Hessian
// [[Rcpp::export]]
mat hessll(vec& b, const colvec& Y, const colvec& lfactY, const mat& X, const mat& Z, const mat& D,
           const rowvec& K, const vec& l0u, const mat& Fu, 
           const vec& beta, const vec& eta, const vec& gr, const int nK){
  mat po = po_hess(b, X, Z, D, beta);
  mat diaggr = diagmat(gr);
  vec kernel = l0u % kron(exp(K * eta), exp(repmat(Fu, 1, nK) * (gr % b)));
  mat diagkern = diagmat(kernel);
  return po + (-diaggr * repmat(Fu, 1, nK).t() * diagkern * repmat(Fu, 1, nK) * diaggr);
}














