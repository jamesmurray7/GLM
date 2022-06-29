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
vec round_N(vec& x, int N){
  // If N is 1000, then round x to hundreths; if X is 10,000 then round x to thousandths.
  double nn = (double)N/10.0;
  return floor(x*nn + 0.50)/nn;
}

//[[Rcpp::export]]
mat gen_lambda_mat(int N, int summax){
  vec mus = linspace(1.0/((double)N/10.0), 10.00, N);
  vec nus = mus;
  mat out = zeros<mat>(N-1, N);
  for(int m = 0; m < (N-1); m ++){
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
  vec mus = linspace(1.0/((double)N/10.0), 10.00, N);
  vec nus = mus;
  mat out = zeros<mat>(N-1, N);
  if((N-1) != lambdamat.n_rows){
    Rcout << "N mistmatch with lambda.n_row" << std::endl;
    return out;
  }else{
    for(int m = 0; m < (N-1); m++){
      for(int n = 0; n < N; n++){
        double lambda = as_scalar(lambdamat(m, n));
        out(m, n) += calc_V(mus[m], lambda, nus[n], summax);
      }
      if (m % 100 == 0) Rcout << m << std::endl;
    }
    return out;
  }
}

// [[Rcpp::export]]
mat gen_logZ_mat(int N, int summax, mat& lambdamat){
  vec mus = linspace(1.0/((double)N/10.0), 10.00, N);
  vec nus = mus;
  mat out = zeros<mat>(N-1, N);
  if((N-1) != lambdamat.n_rows){
    Rcout << "N mistmatch with lambda.n_row" << std::endl;
    return out;
  }else{
    for(int m = 0; m < (N-1); m++){
      for(int n = 0; n < N; n++){
        double log_lambda = log(as_scalar(lambdamat(m, n)));
        out(m, n) += logZ_c_scalar(log_lambda, nus[n], summax);
      }
      if (m % 100 == 0) Rcout << m << std::endl;
    }
    return out;
  }
}

// [[Rcpp::export]]
vec mu_fix(vec &b, int N){
  double small = 1./((double)N/10.0);
  vec x = round_N(b, N);
  return clamp(x, 0.0+small, 10.0-small);
}  

// The joint density ---------------------------------------------------
static double log2pi = log(2.0 * M_PI); 

// [[Rcpp::export]]
double ll_cmp(vec& loglambda, vec& nu, vec& logZ, vec& Y, vec& lY){
  return sum(Y % loglambda - nu % lY - logZ);
}

// [[Rcpp::export]]
vec mat_lookup(vec &m, vec &n, mat& MAT){ // Simply change what matrix is supplied as third argument.
  vec out = vec(m.size());
  for(int i = 0; i < m.size(); i++){
    out[i] += (double)MAT((int)m[i], (int)n[i]);
  }
  return out;
}

// [[Rcpp::export]]
double joint_density(vec& b, mat& X, vec& Y, vec& lY, mat& Z, mat& G,     // Data matrices
                     vec& beta, vec& delta, mat& D,                       // Longit. + RE parameters
                     rowvec& S, mat& SS, rowvec& Fi, mat& Fu, double l0i, rowvec& haz,
                     int Delta, double gamma, vec& zeta, mat& lambdamat, mat& Vmat, mat& logZmat, int N){
  vec mu = exp(X * beta + Z * b);
  vec nu = exp(G * delta);
  vec mu2 = mu_fix(mu, N), nu2 = mu_fix(nu, N);          // Fix \mu s.t. it lies (0.01,10.00). Rounding is controlled by N (1e3 or 1e5).
  // Convert to indices for matrix lookup (+ zero-indexing!!)
  vec m = (mu2/(1./((double)N/10.))) - 1.00, n = (nu2/(1./((double)N/10.))) - 1.00;      
  // Lookup lambda
  vec lambda = mat_lookup(m, n, lambdamat);
  vec loglambda = log(lambda);
  //vec V = getV(m, n, Vmat);
  vec logZ = mat_lookup(m, n, logZmat);
  double ll = ll_cmp(loglambda, nu2, logZ, Y, lY);
  int q = b.size();
  double temp = 0.0;
  if(Delta == 1) temp = log(l0i);
  return -ll + -1.0 * as_scalar(-(double)q/2.0 * log2pi - 0.5 * log(det(D)) - 0.5 * b.t() * D.i() * b + 
                                temp + Delta * (S * zeta + Fi * (gamma * b)) - haz * exp(SS * zeta + Fu * (gamma * b)));
}
  
// [[Rcpp::export]]
vec joint_density_ddb(vec& b, mat& X, vec& Y, vec& lY, mat& Z, mat& G,     // Data matrices
                      vec& beta, vec& delta, mat& D,                       // Longit. + RE parameters
                      rowvec& S, mat& SS, rowvec& Fi, mat& Fu, double l0i, rowvec& haz,
                      int Delta, double gamma, vec& zeta, mat& lambdamat, mat& Vmat, mat& logZmat, int N){
  // Score of CMP
  vec mu = exp(X * beta + Z * b);
  vec nu = exp(G * delta);
  vec mu2 = mu_fix(mu, N), nu2 = mu_fix(nu, N);          // Fix \mu s.t. it lies (0.01,10.00). Rounding is controlled by N (1e3 or 1e5).
  // Convert to indices for matrix lookup (+ zero-indexing!!)
  vec m = (mu2/(1./((double)N/10.))) - 1.00, n = (nu2/(1./((double)N/10.))) - 1.00;    
  vec lambda = mat_lookup(m, n, lambdamat);
  vec loglambda = log(lambda);
  vec V = mat_lookup(m, n, Vmat);
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
                      int Delta, double gamma, vec& zeta, mat& lambdamat, mat& Vmat, mat& logZmat, int N, double eps){
  int q = b.size();
  mat out = zeros<mat>(q, q);
  vec f0 = joint_density_ddb(b, X, Y, lY, Z, G, beta, delta, D, S, SS, Fi, Fu, l0i, haz, Delta, gamma, zeta, lambdamat, Vmat, logZmat, N);
  for(int i = 0; i < q; i++){
    vec bb = b;
    double xi = std::max(1.0, b[i]);
    bb[i] = b[i] + (xi * eps);
    vec fdiff = joint_density_ddb(bb, X, Y, lY, Z, G, beta, delta, D, S, SS, Fi, Fu, l0i, haz, Delta, gamma, zeta, lambdamat, Vmat, logZmat, N) - f0;
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
  
  
/* Alternative update to delta
 * */

double ll_cmp2(vec& loglambda, vec& nu, vec& logZ, vec& Y, vec& lY){
  return sum(Y % loglambda - nu % lY - logZ);
}

// [[Rcpp::export]]
double E_llcmp_delta(vec& delta, mat& G, vec& b, mat& X, mat& Z, vec& Y, vec& lY, vec& beta, vec& tau, vec& w, vec& v, int N,
                     mat& lambdamat, mat& logZmat){
  // nu unchanging
  vec nu = exp(G * delta);
  vec nu2 = mu_fix(nu, N);
  vec eta = X * beta + Z * b;
  // init output
  int gh = w.size();
  double out = 0.0;
  for(int l = 0; l < gh; l++){
    vec this_eta = eta + v[l] * tau;
    vec mu = exp(this_eta);
    vec mu2 = mu_fix(mu, N);
    vec m = (mu2/(1./((double)N/10.))) - 1.00, n = (nu2/(1./((double)N/10.))) - 1.00;    
    vec lambda = mat_lookup(m, n, lambdamat);
    vec loglambda = log(lambda);
    vec logZ = mat_lookup(m, n, logZmat);
    out += w[l] * ll_cmp2(loglambda, nu2, logZ, Y, lY);
  }
  return out;
}

// Forward differencing ----
// [[Rcpp::export]]
vec Sdelta_fdiff(vec& delta, mat& G, vec& b, mat& X, mat& Z, vec& Y, vec& lY, vec& beta, vec& tau, vec& w, vec& v, int N,
                 mat& lambdamat, mat& logZmat, double eps){
  int a = delta.size();
  vec out = vec(a);
  double f0 = E_llcmp_delta(delta, G, b, X, Z, Y, lY, beta, tau, w, v, N, lambdamat, logZmat);
  for(int i = 0; i < a; i++){
    vec aa = delta;
    double xi = std::max(delta[i], 1.0);
    aa[i] = delta[i] + eps * xi;
    double feps = E_llcmp_delta(aa, G, b, X, Z, Y, lY, beta, tau, w, v, N, lambdamat, logZmat) - f0;
    out[i] = feps/(aa[i]-delta[i]);
  }
  return out;
}

// [[Rcpp::export]]
mat Hdelta_fdiff(vec& delta, mat& G, vec& b, mat& X, mat& Z, vec& Y, vec& lY, vec& beta, vec& tau, vec& w, vec& v, int N,
                 mat& lambdamat, mat& logZmat, double eps1, double eps2){
  int a = delta.size();
  mat out = zeros<mat>(a,a);
  vec f0 = Sdelta_fdiff(delta, G, b, X, Z, Y, lY, beta, tau, w, v, N, lambdamat, logZmat, eps1);
  for(int i = 0; i < a; i++){
    vec aa = delta;
    double xi = std::max(delta[i], 1.0);
    aa[i] = delta[i] + eps2 * xi;
    vec feps = Sdelta_fdiff(aa, G, b, X, Z, Y, lY, beta, tau, w, v, N, lambdamat, logZmat, eps1) - f0;
    out.col(i) = feps/(aa[i]-delta[i]);
  }
  return 0.5 * (out + out.t());
}

// Central differencing ----
// [[Rcpp::export]]
vec Sdelta_cdiff(vec& delta, mat& G, vec& b, mat& X, mat& Z, vec& Y, vec& lY, vec& beta, vec& tau, vec& w, vec& v, int N,
                 mat& lambdamat, mat& logZmat, double eps = 1e-2){
  int a = delta.size();
  vec out = vec(a);
  double f0 = E_llcmp_delta(delta, G, b, X, Z, Y, lY, beta, tau, w, v, N, lambdamat, logZmat);
  for(int i = 0; i < a; i++){
    vec aa = delta, bb = delta;
    double xi = std::max(std::fabs(delta[i]), 1.0);
    aa[i] = delta[i] + eps * xi;
    bb[i] = delta[i] - eps * xi;
    double E1 = E_llcmp_delta(aa, G, b, X, Z, Y, lY, beta, tau, w, v, N, lambdamat, logZmat);
    double E2 = E_llcmp_delta(bb, G, b, X, Z, Y, lY, beta, tau, w, v, N, lambdamat, logZmat);
    out[i] = (E1 - E2)/(aa[i]-bb[i]);
  }
  return out;
}