#include <RcppArmadillo.h>
#include <math.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// Helper functions - reshaping vech(A) to A.
// (https://en.wikipedia.org/wiki/Duplication_and_elimination_matrices)
// [[Rcpp::export]]
arma::mat duplication_matrix(const int &n) {
  arma::mat out((n*(n+1))/2, n*n, arma::fill::zeros);
  for (int j = 0; j < n; ++j) {
    for (int i = j; i < n; ++i) {
      arma::vec u((n*(n+1))/2, arma::fill::zeros);
      u(j*n+i-((j+1)*j)/2) = 1.0;
      arma::mat T(n,n, arma::fill::zeros);
      T(i,j) = 1.0;
      T(j,i) = 1.0;
      out += u * arma::trans(arma::vectorise(T));
    }
  }
  return out.t();
}

// [[Rcpp::export]]
mat vech2mat(vec& x, int q){
  vec vecx = duplication_matrix(q) * x;
  return reshape(vecx, q, q);
}

// Functions to calculate Z --------------------------------------------

// [[Rcpp::export]]
arma::vec SEQ_Z(long double summax){ // Shorthand for Z_(lambda_i, nu_i)
  return arma::linspace(0, summax, summax + 1.0);
}

// [[Rcpp::export]]
double calcZ_scalar(double lambda, double nu, int summax){
  double out = 0.0;
  vec js = SEQ_Z(summax);
  for(int j = 0; j < js.size(); j++){
    out += pow(lambda, js[j])/pow(gamma(js[j] + 1.0), nu);
  }
  return out;
}

// [[Rcpp::export]]
vec calcZ(vec& lambda, vec& nu, int summax){
  vec out = vec(lambda.size());
  vec js = SEQ_Z(summax);
  for(int i = 0; i < out.size(); i++){
    for(int j = 0; j < js.size(); j++){
      out[i] += pow(lambda[i], js[j])/pow(gamma(js[j] + 1.0), nu[i]);
    }
  }
  return out;
}

// Function used to find lambda from supplied mu (see _Functions.R) ----

// [[Rcpp::export]]
double mu_lambdaZ_eq(double lambda, double mu, double nu, int summax){
  vec js = SEQ_Z(summax);
  double Z = calcZ_scalar(lambda, nu, summax);
  double rhs = 0.0;
  for(int j = 0; j < js.size(); j++){
    rhs += (js[j] * pow(lambda, js[j])) / (pow(gamma(js[j] + 1.0), nu) * Z);
  }
  return mu - rhs;
}

// Log-likelihood of Y~CMP(lambda, nu) ---------------------------------

// [[Rcpp::export]]
double ll_cmpC(vec& lambda, vec& nu, int summax, vec& Y, vec& lY){
	vec Z = calcZ(lambda, nu, summax);
	double lhs = as_scalar(Y.t() * log(lambda) - nu.t() * lY);
	return lhs - sum(Z);
}

// Log-likelihood of f(b|D) --------------------------------------------

// [[Rcpp::export]]
double ll_b(vec& b, mat& D){
	double q = b.size();
	return as_scalar(
		-q/2.0 * log(2.0 * M_PI) - 0.5 * log(det(D)) - 0.5 * b.t() * D.i() * b
	);
}
