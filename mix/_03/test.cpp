#include <RcppArmadillo.h>
#include <math.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
vec test(vec& v, int start, int end){
  return v.subvec(start, end);
}

/*** R
b <- c(1,2,3,4,5)
test(b,0,1)
test(b,2,2)
test(b,3,4)
*/
