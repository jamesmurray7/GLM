#include <RcppArmadillo.h>
#include <math.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
vec uvectest(vec& x, int start, int end) {
  return x.subvec(start, end);
}


// [[Rcpp::export]]
rowvec newtest(mat& M, int row, uvec& cols){
  rowvec r = M.row(row);
  return r.elem(cols).t();
}

// [[Rcpp::export]]
rowvec submattest(mat& M, uvec row, uvec& cols){
  return M.submat(row, cols);
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
# x <- 1:5
# uvectest(x, 0, 2)
# uvectest(x, 1, 2)

M <- matrix(1:25,5,5,byr=T)
print(M)
newtest(M, 1, c(0, 1))
newtest(M, 4, c(1,2,4))
*/
