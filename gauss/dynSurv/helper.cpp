#include <RcppArmadillo.h>
#include <math.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

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
mat vech2mat(vec& x, const int q){
  vec xx = duplication_matrix(q) * x;
  return reshape(xx, q, q);
}

// [[Rcpp::export]]
long double S_Cpp(List L, vec& gamma, vec& eta, vec& b){
  rowvec l0u = L["l0u.t"];
  mat KK = L["KK.t"];
  mat Fu = L["Fu.t"];
  return as_scalar(exp(-l0u * exp(KK * eta + repmat(Fu, 1, 3) * (gamma % b))));
}

// [[Rcpp::export]]
long double S_Cpp2(List Lt, List Lu, vec& gamma, vec& eta, vec& b){
  rowvec l0u_t = Lt["l0u.t"];
  mat KK_t = Lt["KK.t"];
  mat Fu_t = Lt["Fu.t"];
  rowvec l0u_u = Lu["l0u.t"];
  mat KK_u = Lu["KK.t"];
  mat Fu_u = Lu["Fu.t"];
  return as_scalar(exp(-l0u_u * exp(KK_u * eta + Fu_u * (gamma % b))) / 
                   exp(-l0u_t * exp(KK_t * eta + Fu_t * (gamma % b))));
}



