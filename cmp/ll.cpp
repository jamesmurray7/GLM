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

// credit due to Thomas Fung https://github.com/thomas-fung/mpcmp/blob/main/src/utility.cpp
// A lot easier to simply work on the log scale
// [[Rcpp::export]]
vec logZ_c(vec& log_lambda, vec& nu, int summax) {
  // Control loop
  // int maxiter = 1e4;
  double log_epsilon = std::log(1e-10);
  // Output vector
  int n = log_lambda.size();
  vec out = vec(n);
  // Compute logz
  for (int i = 0; i < n; ++i) {
    double logz  = 0;
    double logz_ = 0;
    for (int j = 1; j < summax; ++j) {
      logz_ += log_lambda[i] - nu[i] * log((double)j);
      logz = R::logspace_add(logz, logz_);
      if (logz_ - logz < log_epsilon) break;
    }
    out[i] = logz;
  }
  return out;
}

// [[Rcpp::export]]
double logZ_c_scalar(double log_lambda, double nu, int summax) {
  // Control loop
  // int maxiter = 1e4;
  double log_epsilon = std::log(1e-10);
  // Output vector
  double out = 0.0;
  // Compute logz
  double logz  = 0;
  double logz_ = 0;
  for (int j = 1; j < summax; ++j) {
    logz_ += log_lambda - nu * log((double)j);
    logz = R::logspace_add(logz, logz_);
    if (logz_ - logz < log_epsilon) break;
  }
  out = logz;
  return out;
}

// double calcZ_scalar(double lambda, double nu, int summax){
//   double out = 0.0;
//   vec js = SEQ_Z(summax);
//   for(int j = 0; j < js.size(); j++){
//     out += pow(lambda, js[j])/pow(tgamma(js[j] + 1.0), nu);
//   }
//   return out;
// }
// 
// vec calcZ(vec& lambda, vec& nu, int summax){
//   vec out = vec(lambda.size());
//   vec js = SEQ_Z(summax);
//   for(int i = 0; i < out.size(); i++){
//     for(int j = 0; j < js.size(); j++){
//       out[i] += pow(lambda[i], js[j])/pow(tgamma(js[j] + 1.0), nu[i]);
//     }
//   }
//   return out;
// }

// Function used to find lambda from supplied mu (see _Functions.R) ----

// [[Rcpp::export]]
double mu_lambdaZ_eq(double lambda, double mu, double nu, int summax){
  vec js = SEQ_Z(summax);
  double Z = exp(logZ_c_scalar(log(lambda), nu, summax));
  double rhs = 0.0;
  for(int j = 0; j < js.size(); j++){
    rhs += (js[j] * pow(lambda, js[j])) / (pow(tgamma(js[j] + 1.0), nu) * Z);
  }
  return mu - rhs;
}

// Log-likelihood of Y~CMP(lambda, nu) ---------------------------------

// [[Rcpp::export]]
double ll_cmpC(vec& lambda, vec& nu, int summax, vec& Y, vec& lY){
  vec loglambda = log(lambda);
	vec Z = exp(logZ_c(loglambda, nu, summax));
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

// And the first derivative wrt b
// [[Rcpp::export]]
vec S_ll_b(vec& b, mat& D){
  return -1.0 * D.i() * b;
}

// Computing Expectations ----------------------------------------------
// [[Rcpp::export]]
vec E_means(vec& lambda, vec& nu, vec& logZ, int summax){
  vec js = SEQ_Z(summax);
  vec out = vec(lambda.size());
  for(int i = 0; i < lambda.size(); i++){
    for(int j = 1; j < js.size(); j++){
      out[i] += exp(log(js[j] - 1.0) + (js[j] - 1.0) *  log(lambda[i]) - nu[i] * lgamma(js[j]) - logZ[i]);
    }
  }
  return out;
}

// [[Rcpp::export]]
vec E_logfactY(vec& lambda, vec& nu, vec& logZ, int summax){
  vec js = SEQ_Z(summax);
  vec out = vec(lambda.size());
  for(int i = 0; i < lambda.size(); i++){
    for(int j = 1; j < js.size(); j++){
      out[i] += lgamma(js[j]) * exp((js[j] - 1) * log(lambda[i]) - nu[i] * lgamma(js[j]) - logZ[i]);
    }
  }
  return out;
}

// [[Rcpp::export]]
vec E_YlogfactY(vec& lambda, vec& nu, vec& logZ, int summax){
  vec js = SEQ_Z(summax);
  vec out = vec(lambda.size());
  for(int i = 0; i < lambda.size(); i++){
    for(int j = 1; j < js.size(); j++){
      out[i] += exp(log(lgamma(js[j])) + log(js[j] - 1.0) + (js[j] - 1.0) * log(lambda[i]) - nu[i] * lgamma(js[j]) - logZ[i]);
    }
  }
  return out;
}

// Variances -----------------------------------------------------------
// [[Rcpp::export]]
vec Var_means(vec& lambda, vec& nu, vec& logZ, int summax){
  vec js = SEQ_Z(summax);
  vec out = vec(lambda.size());
  for(int i = 0; i < lambda.size(); i++){
    for(int j = 1; j < js.size(); j++){
      out[i] += exp(2.0 * log(js[j] - 1.0) + (js[j] - 1.0) * log(lambda[i]) - nu[i] * lgamma(js[j]) - logZ[i]);
    }
  }
  return out - pow(E_means(lambda, nu, logZ, summax), 2.0);
}

// Variance of lfact(Y)
// [[Rcpp::export]]
vec Var_logfactY(vec& lambda, vec& nu, vec& logZ, int summax){
  vec js = SEQ_Z(summax);
  vec out = vec(lambda.size());
  for(int i = 0; i < lambda.size(); i++){
    for(int j = 1; j < js.size(); j++){
        out[i] += pow(lgamma(js[j]), 2.0) * exp((js[j] - 1.0) * log(lambda[i]) - nu[i] * lgamma(js[j]) - logZ[i]);
    }
  }
  return out - pow(E_logfactY(lambda, nu, logZ, summax), 2.0);
}

// Score on \beta ------------------------------------------------------
// Variance V(\mu, Y) 

// [[Rcpp::export]]
double V_scalar(double mu, double lambda, double nu, int summax){
	double out = 0.0;
	double Z = exp(logZ_c_scalar(log(lambda), nu, summax));
	vec js = SEQ_Z(summax);
	for(int j = 0; j < summax; j++){
	  out += pow(js[j] - mu, 2.0) * pow(lambda, js[j]) / (pow(tgamma(js[j] + 1.0), nu) * Z);
	}
	return out;
}

// [[Rcpp::export]]
vec V_mu_lambda(vec& mu, vec& lambda, vec& nu, int summax){
  vec out(lambda.size());
  for(int i = 0; i < lambda.size(); i++){
    out[i] = V_scalar(mu[i], lambda[i], nu[i], summax);
  }
  return out;
}

// A wrapper that uses the above V_scalar and calculates Score(\mu).

// [[Rcpp::export]]
vec Smu(vec& mu, vec& lambda, vec& nu, vec& Y, int summax){
  vec V  = vec(lambda.size());
  for(int v = 0; v < V.size(); v++){
    V[v] = V_scalar(mu[v], lambda[v], nu[v], summax);
  }
  return (Y - mu) / V; 
}





