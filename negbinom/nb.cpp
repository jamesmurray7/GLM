#include <RcppArmadillo.h>
#include <math.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// Negative log-likelihood for joint density f(Y,T|b,...)
// [[Rcpp::export]]
double joint_density(vec& b, mat& X, vec& Y, mat& Z, vec& beta, double theta, mat& D,
                     int Delta, rowvec& K, rowvec Fi, double l0i, 
                     mat& KK, mat& Fu, rowvec& haz, double gamma, vec& eta){
	vec this_eta = X * beta + Z * b;
	
}













