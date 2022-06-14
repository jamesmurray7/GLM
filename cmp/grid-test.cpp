// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <math.h>
#include <float.h>
using namespace Rcpp;
using namespace arma;
#define EPSILON DBL_EPSILON

/*  ------------------------------------------------------------------------ 
 *  1.
 *  Functions required to solve for lambda, the rate parameter via the mean mu.
 *  ------------------------------------------------------------------------ 
*/

arma::vec SEQ_Z(long double summax){ // Shorthand for Z_(lambda_i, nu_i)
  return arma::linspace(0, summax, summax + 1.0);
}

// Normalising constants, Z -----------------------------------------------
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

// [[Rcpp::export]]
double mu_lambdaZ_eq(double lambda, double mu, double nu, int summax){
  vec js = SEQ_Z(summax);
  // double Z = exp(logZ_c_scalar(log(lambda), nu, summax));
  double rhs = 0.0;
  for(int j = 0; j < js.size(); j++){
    rhs += (js[j] - mu) * pow(lambda, js[j]) / pow(tgamma(js[j] + 1.0), nu);
    // rhs += (js[j] * pow(lambda, js[j])) / (pow(tgamma(js[j] + 1.0), nu) * Z);
  }
  return rhs;
}

// UNIROOT ----------------------------------------------------------------
double zeroin_lambda(			/* An estimate of the root */
double ax,				/* Left border | of the range	*/
double bx,				/* Right border| the root is seeked*/
double fa, double fb,		/* f(a), f(b) */
double mu,				    /* Add'l info passed on to f	*/
double nu,            /*    ------ "" ----- */
int summax,           /*    ------ "" ----- */
double Tol,			/* Acceptable tolerance		*/
int Maxit)				/* Max # of iterations */
{
  double a,b,c, fc;			/* Abscissae, descr. see above,  f(c) */
  double tol;
  int maxit;
  
  a = ax;  b = bx;
  c = a;   fc = fa;
  maxit = Maxit + 1; tol =  Tol;
  
  /* First test if we have found a root at an endpoint */
  if(fa == 0.0) {
    Tol = 0.0;
    Maxit = 0;
    return a;
  }
  if(fb ==  0.0) {
    Tol = 0.0;
    Maxit = 0;
    return b;
  }
  
  while(maxit--)		/* Main iteration loop	*/
  {
    double prev_step = b-a;		/* Distance from the last but one
   to the last approximation	*/
  double tol_act;			/* Actual tolerance		*/
  double p;			/* Interpolation step is calcu- */
  double q;			/* lated in the form p/q; divi-
   * sion operations is delayed
   * until the last moment	*/
  double new_step;		/* Step at this iteration	*/
  
  if( fabs(fc) < fabs(fb) )
  {				/* Swap data for b to be the	*/
  a = b;  b = c;  c = a;	/* best approximation		*/
  fa=fb;  fb=fc;  fc=fa;
  }
  tol_act = 2*EPSILON*fabs(b) + tol/2;
  new_step = (c-b)/2;
  
  if( fabs(new_step) <= tol_act || fb == (double)0 )
  {
    Maxit -= maxit;
    Tol = fabs(c-b);
    return b;			/* Acceptable approx. is found	*/
  }
  
  /* Decide if the interpolation can be tried	*/
  if( fabs(prev_step) >= tol_act	/* If prev_step was large enough*/
  && fabs(fa) > fabs(fb) ) {	/* and was in true direction,
   * Interpolation may be tried	*/
  register double t1,cb,t2;
    cb = c-b;
    if( a==c ) {		/* If we have only two distinct	*/
  /* points linear interpolation	*/
  t1 = fb/fa;		/* can only be applied		*/
  p = cb*t1;
  q = 1.0 - t1;
    }
    else {			/* Quadric inverse interpolation*/
  
  q = fa/fc;  t1 = fb/fc;	 t2 = fb/fa;
  p = t2 * ( cb*q*(q-t1) - (b-a)*(t1-1.0) );
  q = (q-1.0) * (t1-1.0) * (t2-1.0);
    }
    if( p>(double)0 )		/* p was calculated with the */
  q = -q;			/* opposite sign; make p positive */
  else			/* and assign possible minus to	*/
  p = -p;			/* q				*/
  
  if( p < (0.75*cb*q-fabs(tol_act*q)/2) /* If b+p/q falls in [b,c]*/
  && p < fabs(prev_step*q/2) )	/* and isn't too large	*/
  new_step = p/q;			/* it is accepted
   * If p/q is too large then the
   * bisection procedure can
   * reduce [b,c] range to more
   * extent */
  }
  
  if( fabs(new_step) < tol_act) {	/* Adjust the step to be not less*/
  if( new_step > (double)0 )	/* than tolerance		*/
  new_step = tol_act;
  else
    new_step = -tol_act;
  }
  a = b;	fa = fb;			/* Save the previous approx. */
  b += new_step;	fb = mu_lambdaZ_eq(b, mu, nu, summax);	/* Do step to a new approxim. */
  if( (fb > 0 && fc > 0) || (fb < 0 && fc < 0) ) {
    /* Adjust c for it to have a sign opposite to that of b */
    c = a;  fc = fa;
  }
  
  }
  /* failed! */
  Tol = fabs(c-b);
  Maxit = -1;
  return b;
}

//[[Rcpp::export]]
double lambda_uniroot(double lower, double upper, double mu, double nu, int summax){
  double f_lower = mu_lambdaZ_eq(lower, mu, nu, summax);
  double f_upper = mu_lambdaZ_eq(upper, mu, nu, summax);
  double lambda = zeroin_lambda(lower, upper, f_lower, f_upper, mu, nu, summax, 1e-4, 100);
  return lambda;
}

//[[Rcpp::export]]
vec lambda_uniroot_wrap(double lower, double upper, vec& mu, vec& nu, int summax){
  int m = mu.size();
  vec lambda = vec(m);
  for(int i = 0; i < m; i++){
    double f_lower = mu_lambdaZ_eq(lower, mu[i], nu[i], summax);
    double f_upper = mu_lambdaZ_eq(upper, mu[i], nu[i], summax);
    lambda[i] = zeroin_lambda(lower, upper, f_lower, f_upper, mu[i], nu[i], summax, 1e-4, 100);
  }
  return lambda;
}

// [[Rcpp::export]]
double calc_V(double mu, double lambda, double nu, int summax){
  double out = 0.0;
  double Z = exp(logZ_c_scalar(log(lambda), nu, summax));
  vec js = SEQ_Z(summax);
  for(int j = 0; j < summax; j++){
    out += pow(js[j] - mu, 2.0) * pow(lambda, js[j]) / (pow(tgamma(js[j] + 1.0), nu) * Z);
  }
  return out;
}

//[[Rcpp::export]]
vec round_hundreth(vec& x){
  return floor(x*100.0 + 0.50)/100.0;
}

//[[Rcpp::export]]
mat gen_lambda_mat(int N, int summax){
  vec mus = linspace(0.01, 10.00, N);
  vec nus = mus;
  mat out = zeros<mat>(N, N);
  for(int m = 0; m < N; m ++){
    for(int n = 0; n < N; n++){
      double l = lambda_uniroot(1e-6, 1e3, mus[m], nus[n], summax);
      if(l == 1e-6){
        out(m, n) += out(m-1, n);
      }else{
        out(m, n) += l;  
      }
    }
    if (m % 100 == 0) Rcout << m << std::endl;
  }
  return out;
}

//[[Rcpp::export]]
mat gen_V_mat(int N, int summax, mat& lambdamat){
  vec mus = linspace(0.01, 10.00, N);
  vec nus = mus;
  mat out = zeros<mat>(N, N);
  if(N != lambdamat.n_rows){
    Rcout << "N mistmatch with lambda.n_row" << std::endl;
    return out;
  }else{
    for(int m = 0; m < N; m++){
      for(int n = 0; n < N; n++){
        out(m, n) += calc_V(mus[m], lambdamat(m, n), nus[n], summax);
      }
      if (m % 100 == 0) Rcout << m << std::endl;
    }
    return out;
  }
}

// [[Rcpp::export]]
mat gen_logZ_mat(int N, int summax, mat& lambdamat){
  vec mus = linspace(0.01, 10.00, N);
  vec nus = mus;
  mat out = zeros<mat>(N, N);
  if(N != lambdamat.n_rows){
    Rcout << "N mistmatch with lambda.n_row" << std::endl;
    return out;
  }else{
    for(int m = 0; m < N; m++){
      for(int n = 0; n < N; n++){
        out(m, n) += logZ_c_scalar(log(lambdamat(m, n)), nus[n], summax);
      }
      if (m % 100 == 0) Rcout << m << std::endl;
    }
    return out;
  }
}

// [[Rcpp::export]]
vec mu_fix(vec &b){
  vec x = round_hundreth(b);
  return clamp(x, 0.010, 10.000);
}  

// The joint density ---------------------------------------------------
static double log2pi = log(2.0 * M_PI); 

// [[Rcpp::export]]
double ll_cmp(vec& loglambda, vec& nu, vec& logZ, vec& Y, vec& lY){
  return as_scalar(Y.t() * loglambda - nu.t() * lY) - sum(logZ);
}

// [[Rcpp::export]]
vec getlambda(vec &m, vec &n, mat& lambdamat){
  vec out = vec(m.size());
  for(int i = 0; i < m.size(); i++){
    out[i] = lambdamat((int)m[i], (int)n[i]);
  }
  return out;
}

// [[Rcpp::export]]
vec getV(vec &m, vec &n, mat& Vmat){
  vec out = vec(m.size());
  for(int i = 0; i < m.size(); i++){
    out[i] = Vmat((int)m[i], (int)n[i]);
  }
  return out;
}

// [[Rcpp::export]]
vec getlogZ(vec &m, vec &n, mat& logZmat){
  vec out = vec(m.size());
  for(int i = 0; i < m.size(); i++){
    out[i] = logZmat((int)m[i], (int)n[i]);
  }
  return out;
}

// [[Rcpp::export]]
double joint_density(vec& b, mat& X, vec& Y, vec& lY, mat& Z, mat& G,     // Data matrices
                     vec& beta, vec& delta, mat& D,                       // Longit. + RE parameters
                     rowvec& S, mat& SS, rowvec& Fi, mat& Fu, double l0i, rowvec& haz,
                     int Delta, double gamma, vec& zeta, mat& lambdamat, mat& Vmat, mat& logZmat){
  vec mu = exp(X * beta + Z * b);
  vec nu = exp(G * delta);
  vec mu2 = mu_fix(mu);                            // Fix \mu s.t. it lies [0.01,10.00]. Shouldn't be necessary with \nu
  vec m = (mu2/0.01) - 0.01, n = (nu/0.01) - 0.01; // Convert to indices for matrix lookup (+ zero-indexing!!)
  vec lambda = getlambda(m, n, lambdamat);
  vec loglambda = log(lambda);
  //vec V = getV(m, n, Vmat);
  vec logZ = getlogZ(m, n, logZmat);
  double ll = ll_cmp(loglambda, nu, logZ, Y, lY);
  int q = b.size();
  double temp = 0.0;
  if(Delta == 1) temp = log(l0i);
  return -ll + -1.0 * as_scalar(-q/2.0 * log2pi - 0.5 * log(det(D)) - 0.5 * b.t() * D.i() * b + 
                                temp + Delta * (S * zeta + Fi * (gamma * b)) - haz * exp(SS * zeta + Fu * (gamma * b)));
}
  
// [[Rcpp::export]]
vec joint_density_ddb(vec& b, mat& X, vec& Y, vec& lY, mat& Z, mat& G,     // Data matrices
                      vec& beta, vec& delta, mat& D,                       // Longit. + RE parameters
                      rowvec& S, mat& SS, rowvec& Fi, mat& Fu, double l0i, rowvec& haz,
                      int Delta, double gamma, vec& zeta, mat& lambdamat, mat& Vmat, mat& logZmat){
  // Score of CMP
  vec mu = exp(X * beta + Z * b);
  vec nu = exp(G * delta);
  vec mu2 = mu_fix(mu);                            // Fix \mu s.t. it lies [0.01,10.00]. Shouldn't be necessary with \nu
  vec m = (mu2/0.01) - 0.01, n = (nu/0.01) - 0.01; // Convert to indices for matrix lookup (+ zero-indexing!!)
  vec lambda = getlambda(m, n, lambdamat);
  vec loglambda = log(lambda);
  vec V = getV(m, n, Vmat);
  //vec logZ = getlogZ(m, n, logZmat);
  mat lhs_mat = diagmat(mu2) * Z;            // Could just be Z.t() * y-mu/V, but think it should be this.
  vec Score = lhs_mat.t() * ((Y-mu2) / V);   //                  **                                     
  
  return -Score + -1.0 * (-D.i() * b + Delta * (Fi.t() * gamma) - 
                          gamma * (Fu.t() * (haz.t() % exp(SS * zeta + Fu * (gamma * b)))));
    
}  

// [[Rcpp::export]]
mat joint_density_sdb(vec& b, mat& X, vec& Y, vec& lY, mat& Z, mat& G,     // Data matrices
                      vec& beta, vec& delta, mat& D,                       // Longit. + RE parameters
                      rowvec& S, mat& SS, rowvec& Fi, mat& Fu, double l0i, rowvec& haz,
                      int Delta, double gamma, vec& zeta, mat& lambdamat, mat& Vmat, mat& logZmat, double eps){
  int q = b.size();
  mat out = zeros<mat>(q, q);
  vec f0 = joint_density_ddb(b, X, Y, lY, Z, G, beta, delta, D, S, SS, Fi, Fu, l0i, haz, Delta, gamma, zeta, lambdamat, Vmat, logZmat);
  for(int i = 0; i < q; i++){
    vec bb = b;
    double xi = std::max(1.0, b[i]);
    bb[i] = b[i] + (xi * eps);
    vec fdiff = joint_density_ddb(bb, X, Y, lY, Z, G, beta, delta, D, S, SS, Fi, Fu, l0i, haz, Delta, gamma, zeta, lambdamat, Vmat, logZmat) - f0;
    out.col(i) = fdiff/(bb[i]-b[i]);
  }
  return 0.5 * (out + out.t()); // Ensure symmetry
}

// Functions associated with the update to \delta ----------------------
// Rewrite ABC calculations in C++ if needed (i.e. it is a bottleneck)

// [[Rcpp::export]]
mat getW2(List ABC, vec& V, vec& nu, mat& G){
  int a = G.n_cols, m = G.n_rows;
  vec A = ABC["A"];
  vec B = ABC["B"];
  vec C = ABC["C"];
  mat out = zeros<mat>(a, a);
  for(int g = 0; g < m; g++){
    out += ((-pow(A[g], 2.0)/V[g]+C[g]) * pow(nu[g], 2.0)) * G.row(g).t() * G.row(g);
  }
  return -1.0 * out;
}

// Defining updates for the survival pair (gamma, zeta) ----------------
// Define the conditional expectation and then take Score AND Hessian via forward differencing
double Egammazeta(vec& gammazeta, vec& b, mat& Sigma,
                  rowvec& S, mat& SS, mat& Fu, rowvec& Fi, vec& haz, int Delta, vec& w, vec& v){
  double g = as_scalar(gammazeta.head(1)); // gamma will always be scalar and the first element
  vec z = gammazeta.subvec(1, gammazeta.size() - 1);  // with the rest of the vector constructed by zeta
  // determine tau
  vec tau = pow(g, 2.0) * diagvec(Fu * Sigma * Fu.t());
  double rhs = 0.0;
  for(int l = 0; l < w.size(); l++){
    rhs += w[l] * as_scalar(haz.t() * exp(SS * z + Fu * (g * b) + v[l] * pow(tau, 0.5)));
  }
  return as_scalar(Delta * (S * z + Fi * (g * b)) - rhs);
}
// [[Rcpp::export]]
vec Sgammazeta(vec& gammazeta, vec& b, mat& Sigma,
               rowvec& S, mat& SS, mat& Fu, rowvec& Fi, vec& haz, int Delta, vec& w, vec& v, long double eps){
  vec out = vec(gammazeta.size());
  double f0 = Egammazeta(gammazeta, b, Sigma, S, SS, Fu, Fi, haz, Delta, w, v);
  for(int i = 0; i < gammazeta.size(); i++){
    vec ge = gammazeta;
    double xi = std::max(ge[i], 1.0);
    ge[i] = gammazeta[i] + xi * eps;
    double fdiff = Egammazeta(ge, b, Sigma, S, SS, Fu, Fi, haz, Delta, w, v) - f0;
    out[i] = fdiff/(ge[i]-gammazeta[i]);
  }
  return out;
}
// [[Rcpp::export]]
mat Hgammazeta(vec& gammazeta, vec& b, mat& Sigma,
               rowvec& S, mat& SS, mat& Fu, rowvec& Fi, vec& haz, int Delta, vec& w, vec& v, long double eps){
  mat out = zeros<mat>(gammazeta.size(), gammazeta.size());
  vec f0 = Sgammazeta(gammazeta, b, Sigma, S, SS, Fu, Fi, haz, Delta, w, v, eps);
  for(int i = 0; i < gammazeta.size(); i++){
    vec ge = gammazeta;
    double xi = std::max(ge[i], 1.0);
    ge[i] = gammazeta[i] + xi * eps;
    vec fdiff = Sgammazeta(ge, b, Sigma, S, SS, Fu, Fi, haz, Delta, w, v, eps) - f0;
    out.col(i) = fdiff/(ge[i]-gammazeta[i]);
  }
  return 0.5 * (out + out.t());
}

// Defining update to the baseline hazard lambda_0 ---------------------
// [[Rcpp::export]]
mat lambdaUpdate(List survtimes, mat& ft, double gamma, vec& zeta,
                 List S, List Sigma, List b, vec& w, vec& v){
  int gh = w.size();
  int id = b.size();
  mat store = zeros<mat>(ft.n_rows, id);
  for(int i = 0; i < id; i++){
    vec survtimes_i = survtimes[i];
    mat Sigma_i = Sigma[i];
    vec b_i = b[i];
    rowvec S_i = S[i];
    for(int j = 0; j < survtimes_i.size(); j++){
      rowvec Fst = ft.row(j);
      double tau = as_scalar(pow(gamma, 2.0) * Fst * Sigma_i * Fst.t());
      vec rhs = gamma * b_i; //vec(b_i.size());
      double mu = as_scalar(exp(S_i * zeta + Fst * rhs));
      for(int l = 0; l < gh; l++){
        store(j, i) += as_scalar(w[l] * mu * exp(v[l] * sqrt(tau)));
      }
    }
  }
  return store;
}
  