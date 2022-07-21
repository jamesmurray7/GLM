// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <math.h>
#include <float.h>
using namespace Rcpp;
using namespace arma;
#define EPSILON DBL_EPSILON

arma::vec SEQ_Z(long double summax){ // Shorthand for Z_(lambda_i, nu_i)
  return arma::linspace(0, summax, summax + 1.0);
}


// Normalising constants, Z -----------------------------------------------
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
double logZ_c_scalar(double log_lambda, double nu, int summax) { // SCALAR VERSION.
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

// Specifically for obtaining pmf in simulations --------------------------
// [[Rcpp::export]]
vec cmp_pmf_scalar(vec& Y, double lambda, double nu, int summax){
  vec out = vec(Y.size());
  double logZ = logZ_c_scalar(log(lambda), nu, summax);
  vec L = exp(Y * log(lambda) - nu * lgamma(Y + 1.0) - logZ);
  L.replace(datum::inf, 1e100);
  return L;
}

// The variance, V --------------------------------------------------------
// [[Rcpp::export]]
double calc_V(double mu, double lambda, double nu, double logZ, int summax){
  double out = 0.0;
  double Z = trunc_exp(logZ);
  vec js = linspace(1, summax, summax);
  for(int j = 0; j < summax; j++){
    out += pow(js[j] - mu, 2.0) * pow(lambda, js[j]) / (pow(tgamma(js[j] + 1.0), nu) * Z);
  }
  return out;
}

vec calc_V_vec(vec& mu, vec& lambda, vec& nu, vec& logZ, int summax){
  int n = mu.size();
  vec out = vec(n);
  for(int i = 0; i < n; i++){
    out[i] = calc_V(mu[i], lambda[i], nu[i], logZ[i], summax);
  }
  return out;
}

// UNIROOT ----------------------------------------------------------------
double mu_lambdaZ_eq(double lambda, double mu, double nu, int summax){
  vec js = SEQ_Z(summax);
  // double Z = exp(logZ_c_scalar(log(lambda), nu, summax));
  double rhs = 0.0;
  for(int j = 0; j < js.size(); j++){
    rhs += (js[j] - mu) * pow(lambda, js[j]) / pow(tgamma(js[j] + 1.0), nu);
  }
  return rhs;
}

double zeroin_lambda(			/* An estimate of the root */
double ax,				        /* Left border | of the range	*/
double bx,				        /* Right border| the root is seeked*/
double fa, double fb,		  /* f(a), f(b) */
double mu,				        /* Add'l info passed on to f	*/
double nu,                /*    ------ "" ----- */
int summax,               /*    ------ "" ----- */
double Tol,			          /* Acceptable tolerance		*/
int Maxit)				        /* Max # of iterations */
{
  double a,b,c, fc;			  /* Abscissae, descr. see above,  f(c) */
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

// Creating lookup matrices -----------------------------------------------
//[[Rcpp::export]]
double lambda_uniroot(double lower, double upper, double mu, double nu, int summax){
  double f_lower = mu_lambdaZ_eq(lower, mu, nu, summax);
  double f_upper = mu_lambdaZ_eq(upper, mu, nu, summax);
  double lambda = zeroin_lambda(lower, upper, f_lower, f_upper, mu, nu, summax, 1e-4, 100);
  return lambda;
}

// [[Rcpp::export]]
mat gen_lambda_mat(int N, int summax){
  vec mus = linspace(1.0/((double)N/10.0), 10.00, N), nus = mus;
  mat out = zeros<mat>(N-1, N);
  for(int m = 0; m < (N-1); m ++){
    for(int n = 0; n < N; n++){
      double init = pow(mus[m] + (nus[n] - 1.0) / (2.0 * nus[n]), nus[n]);
      double lower = std::max(1e-6, init - 10.0);
      double upper = init + 10.0;
      double l = lambda_uniroot(lower, upper, mus[m], nus[n], summax);
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

//[[Rcpp::export]]
mat gen_V_mat(int N, int summax, mat& lambdamat, mat& logZmat){
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
        double logZ = as_scalar(logZmat(m, n));
        out(m, n) += calc_V(mus[m], lambda, nus[n], logZ, summax);
      }
      if (m % 100 == 0) Rcout << m << std::endl;
    }
    return out;
  }
}


// [[Rcpp::export]]
vec lambda_appx(vec& mu, vec& nu, int summax){
  int m = mu.size();
  vec lambda = vec(m);
  vec inits = trunc_exp(nu % log(mu + (nu-1.0)/(2.0 * nu)));
  for(int i = 0; i < m; i++){
    double lower = std::max(inits[i] - 10.0, 1e-6);
    double upper = inits[i] + 10.0;
    double f_lower = mu_lambdaZ_eq(lower, mu[i], nu[i], summax);
    double f_upper = mu_lambdaZ_eq(upper, mu[i], nu[i], summax);
    lambda[i] = zeroin_lambda(lower, upper, f_lower, f_upper, mu[i], nu[i], summax, 1e-4, 100);
  }
  return lambda;
}
  
vec round_N(vec& x, int N){
  // If N is 1000, then round x to hundreths; if X is 10,000 then round x to thousandths.
  double nn = (double)N/10.0;
  return floor(x*nn + 0.50)/nn;
}

vec mu_fix(vec &b, int N){
  double small = 1./((double)N/10.0);
  vec x = round_N(b, N);
  return clamp(x, 0.0+small, 10.0-small);
}  

// Create vectors for lambda, logZ and V given {mu, nu} -------------------
// [[Rcpp::export]]
vec mat_lookup(vec &m, vec &n, mat& MAT){ // Simply change what matrix is supplied as third argument.
  vec out = vec(m.size());
  for(int i = 0; i < m.size(); i++){
    out[i] += (double)MAT((int)m[i], (int)n[i]);
  }
  return out;
}

// [[Rcpp::export]]
vec get_lambda(vec& mu, vec& nu, int summax, int N, mat& lambda_mat){
  
  // Define uvecs for lookups >/< value of mu=10.0
  uvec ge10 = find(mu >= 10.0), lt10 = find(mu < 10.0);
  int a = mu.size(), b = ge10.size(), c = lt10.size();
  vec out = vec(a);
  
  // If there are elements >= 10, fix using the approximaton for lambda...
  if(b > 0){
    vec mu_temp = mu.elem(ge10), nu_temp = nu.elem(ge10);
    Rcout << "ge10" << mu_temp << std::endl;
    out.elem(ge10) = lambda_appx(mu_temp, nu_temp, summax);
  }
  
  // Simply take matrix lookups for remaining elements
  if(c > 0){
    vec mu_temp = mu.elem(lt10), nu_temp = nu.elem(lt10);
    vec mu2 = mu_fix(mu_temp, N), nu2 = mu_fix(nu_temp, N);         
    vec m = (mu2/(1./((double)N/10.))) - 1.00, n = (nu2/(1./((double)N/10.))) - 1.00;      
    out.elem(lt10) = mat_lookup(m, n, lambda_mat);
  }
  
  return out;
}

// [[Rcpp::export]]
vec get_logZ(vec& mu, vec& nu, int summax, int N, vec& lambda, mat& logZ_mat){
  
  // Define lookups for >/< values of mu = 10.0
  uvec ge10 = find(mu >= 10.0), lt10 = find(mu < 10.0);
  int a = mu.size(), b = ge10.size(), c = lt10.size();
  vec out = vec(a);
  
  // If there's elements >= 10, find using the approximation from \lambda.
  if(b > 0){
    vec mu_temp = mu.elem(ge10), nu_temp = nu.elem(ge10), lam_temp = lambda.elem(ge10);
    vec ll = log(lam_temp);
    vec logZ = logZ_c(ll, nu_temp, summax);
    out.elem(ge10) += logZ;
  }
  
  // And simply lookup logZ matrix for remaining mus
  if(c > 0){
    vec mu_temp = mu.elem(lt10), nu_temp = nu.elem(lt10);
    vec mu2 = mu_fix(mu_temp, N), nu2 = mu_fix(nu_temp, N);         
    vec m = (mu2/(1./((double)N/10.))) - 1.00, n = (nu2/(1./((double)N/10.))) - 1.00;      
    out.elem(lt10) = mat_lookup(m, n, logZ_mat);
  }
  
  return out;
}

// [[Rcpp::export]]
vec get_V(vec& mu, vec& nu, int summax, int N, vec& lambda, vec& logZ, mat& V_mat){
  
  // Define lookups for >/< values of mu = 10.0
  uvec ge10 = find(mu >= 10.0), lt10 = find(mu < 10.0);
  int a = mu.size(), b = ge10.size(), c = lt10.size();
  vec out = vec(a);
  
  // If there's elements >= 10, find using the value from logZ and lambda.
  if(b > 0){
    vec mu_temp = mu.elem(ge10), nu_temp = nu.elem(ge10), lZ_temp = logZ.elem(ge10), lam_temp = lambda.elem(ge10);
    out.elem(ge10) = calc_V_vec(mu_temp, lam_temp, nu_temp, lZ_temp, summax);
  }
  
  // And simply lookup V matrix for remaining mus
  if(c > 0){
    vec mu_temp = mu.elem(lt10), nu_temp = nu.elem(lt10);
    vec mu2 = mu_fix(mu_temp, N), nu2 = mu_fix(nu_temp, N);         
    vec m = (mu2/(1./((double)N/10.))) - 1.00, n = (nu2/(1./((double)N/10.))) - 1.00;      
    out.elem(lt10) = mat_lookup(m, n, V_mat);
  }
  
  return out;
}

// Defining a joint density -----------------------------------------------
static double log2pi = log(2.0 * M_PI); 

// [[Rcpp::export]]
double ll_cmp(vec& loglambda, vec& nu, vec& logZ, vec& Y, vec& lY){
  return sum(Y % loglambda - nu % lY - logZ);
}

// [[Rcpp::export]]
double joint_density(vec& b, mat& X, vec& Y, vec& lY, mat& Z, mat& G,     // Data matrices
                     vec& beta, vec& delta, mat& D,                       // Longit. + RE parameters
                     rowvec& S, mat& SS, rowvec& Fi, mat& Fu, double l0i, rowvec& haz,
                     int Delta, double gamma, vec& zeta, mat& lambdamat, mat& Vmat, mat& logZmat, int N, int summax){
  // Define mu, nu
  vec mu = exp(X * beta + Z * b);
  vec nu = exp(G * delta);
  // Lookup lambda and logZ
  // vec lambda = get_lambda(mu, nu, summax, N, lambdamat);
  // vec loglambda = log(lambda);
  // vec logZ = get_logZ(mu, nu, summax, N, lambda, logZmat);
  vec lambda = lambda_appx(mu, nu, summax);
  vec loglambda = log(lambda);
  vec logZ = logZ_c(loglambda, nu, summax);
  // Calculate loglik CMP and then for other consituent distns.
  double ll = ll_cmp(loglambda, nu, logZ, Y, lY);
  int q = b.size();
  double temp = 0.0;
  if(Delta == 1) temp = log(l0i);
  return -ll + -1.0 * as_scalar(-(double)q/2.0 * log2pi - 0.5 * log(det(D)) - 0.5 * b.t() * D.i() * b + 
                                temp + Delta * (S * zeta + Fi * (gamma * b)) - haz * exp(SS * zeta + Fu * (gamma * b)));
}

// Its first derivative wrt b ----
// [[Rcpp::export]]
vec joint_density_ddb(vec& b, mat& X, vec& Y, vec& lY, mat& Z, mat& G,     // Data matrices
                      vec& beta, vec& delta, mat& D,                       // Longit. + RE parameters
                      rowvec& S, mat& SS, rowvec& Fi, mat& Fu, double l0i, rowvec& haz,
                      int Delta, double gamma, vec& zeta, mat& lambdamat, mat& Vmat, mat& logZmat, int N, int summax){
  // Define mu, nu
  vec mu = exp(X * beta + Z * b);
  // Rcout << "mu: " << mu << std::endl;
  vec nu = exp(G * delta);
  // Lookup lambda, logZ and V
  // vec lambda = get_lambda(mu, nu, summax, N, lambdamat);
  // vec logZ = get_logZ(mu, nu, summax, N, lambda, logZmat); // Need this to work out V.
  // vec V = get_V(mu, nu, summax, N, lambda, logZ, Vmat);
  vec lambda = lambda_appx(mu, nu, summax);
  vec loglambda = log(lambda);
  vec logZ = logZ_c(loglambda, nu, summax);
  vec V = calc_V_vec(mu, lambda, nu, logZ, summax);
  // Rcout << "V: " << V << std::endl;
  // Rcout << "lambda: " << lambda << std::endl;
  // Rcout << "logZ: " << logZ << std::endl;
  
  // Score of CMP and other constituent distns.
  mat lhs_mat = diagmat(mu) * Z;          
  vec Score = lhs_mat.t() * ((Y-mu) / V); 
  // Rcout << "CMP Score: " << Score << std::endl;
  
  return -Score + -1.0 * (-D.i() * b + Delta * (Fi.t() * gamma) - 
                          gamma * (Fu.t() * (haz.t() % exp(SS * zeta + Fu * (gamma * b)))));
  
}  

// And take its second derivative via forward differencing ----
// [[Rcpp::export]]
mat joint_density_sdb(vec& b, mat& X, vec& Y, vec& lY, mat& Z, mat& G,     // Data matrices
                      vec& beta, vec& delta, mat& D,                       // Longit. + RE parameters
                      rowvec& S, mat& SS, rowvec& Fi, mat& Fu, double l0i, rowvec& haz,
                      int Delta, double gamma, vec& zeta, mat& lambdamat, mat& Vmat, mat& logZmat, int N, int summax, double eps){
  int q = b.size();
  mat out = zeros<mat>(q, q);
  // Rcout << "0------------" << std::endl;
  vec f0 = joint_density_ddb(b, X, Y, lY, Z, G, beta, delta, D, S, SS, Fi, Fu, l0i, haz, Delta, gamma, zeta, lambdamat, Vmat, logZmat, N, summax);
  for(int i = 0; i < q; i++){
    // Rcout << q << "-------------" << std::endl;
    vec bb = b;
    double xi = std::max(1.0, b[i]);
    bb[i] = b[i] + (xi * eps);
    vec fdiff = joint_density_ddb(bb, X, Y, lY, Z, G, beta, delta, D, S, SS, Fi, Fu, l0i, haz, Delta, gamma, zeta, lambdamat, Vmat, logZmat, N, summax) - f0;
    out.col(i) = fdiff/(bb[i]-b[i]);
  }
  return 0.5 * (out + out.t()); // Ensure symmetry
}

// Updates for the survival pair (gamma, zeta) ----------------------------

// Define the conditional expectation
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

// Then take Score AND Hessian via forward differencing.
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

// Defining update to the baseline hazard lambda_0 ------------------------
// [[Rcpp::export]]
mat lambdaUpdate(List survtimes, mat& ft, double gamma, vec& zeta,
                 List S, List Sigma, List b, vec& w, vec& v){
  int gh = w.size(), id = b.size();
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

// Creation of elements A, B, C -------------------------------------------
vec ElY(vec& lambda, vec& nu, vec& Z, int summax){
  int l = lambda.size();
  mat M = zeros<mat>(l, summax + 1);
  for(int j = 1; j <= summax; j++){
    M.col(j) += lgamma((double)j) * exp(((double)j-1.0) * log(lambda) - nu * lgamma((double)j) - Z);
  }
  return sum(M, 1); // Rowsums!
}

vec EYlY(vec& lambda, vec& nu, vec& Z, int summax){
  int l = lambda.size();
  mat M = zeros<mat>(l, summax + 1);
  for(int j = 1; j <= summax; j++){
    M.col(j) += exp(
      log((double)j - 1.0) + log(lgamma((double)j)) + ((double)j - 1.0) * log(lambda) - nu * lgamma((double)j) - Z
    );
  }
  return sum(M, 1); // Rowsums!
}

vec VlY(vec& lambda, vec& nu, vec& Z, vec& B, int summax){
  int l = lambda.size();
  mat M = zeros<mat>(l, summax + 1);
  for(int j = 1; j <= summax; j++){
    M.col(j) += pow(lgamma((double)j), 2.0) * exp( ((double)j - 1.0) * log(lambda) - nu * lgamma((double)j) - Z );
  }
  return sum(M, 1) - pow(B, 2.0); // Rowsums!
}

// [[Rcpp::export]]
List A_B_C(vec& b, mat& X, mat& Z, vec& beta, vec& delta, mat& G, vec & lambda, vec& logZ, int summax, int N){
  vec mu = exp(X * beta + Z * b);
  vec nu = exp(G * delta);
  // Calculate B, C and A 
  vec B = ElY(lambda, nu, logZ, summax);
  vec C = VlY(lambda, nu, logZ, B, summax);
  vec A = EYlY(lambda, nu, logZ, summax) - mu % B;
  List out = List::create(Named("A") = A , _["B"] = B, _["C"] = C);
  return out;
}

// Update to the dispersion parameter \delta ------------------------------
// The conditional expectation on \b.
// [[Rcpp::export]]
double E_llcmp_delta(vec& delta, mat& G, vec& b, mat& X, mat& Z, vec& Y, vec& lY, vec& beta, vec& tau, vec& w, vec& v, int summax){
  // nu unchanging
  vec nu = exp(G * delta);
  vec eta = X * beta + Z * b;
  // init output
  int gh = w.size();
  double out = 0.0;
  // Loop over quad. nodes.
  for(int l = 0; l < gh; l++){
    vec this_eta = eta + v[l] * tau;
    vec mu = exp(this_eta);
    vec lambda = lambda_appx(mu, nu, summax);
    vec loglambda = log(lambda);
    vec logZ = logZ_c(loglambda, nu, summax);
    out += w[l] * ll_cmp(loglambda, nu, logZ, Y, lY);
  }
  return out;
}

// The score via central differencing
// [[Rcpp::export]]
vec Sdelta_cdiff(vec& delta, mat& G, vec& b, mat& X, mat& Z, vec& Y, vec& lY, vec& beta, vec& tau, vec& w, vec& v, int summax, double eps){
  int a = delta.size();
  vec out = vec(a);
  for(int i = 0; i < a; i++){
    vec aa = delta, bb = delta;
    double xi = std::max(std::fabs(delta[i]), 1.0);
    aa[i] = delta[i] + eps * xi;
    bb[i] = delta[i] - eps * xi;
    double E1 = E_llcmp_delta(aa, G, b, X, Z, Y, lY, beta, tau, w, v, summax);
    double E2 = E_llcmp_delta(bb, G, b, X, Z, Y, lY, beta, tau, w, v, summax);
    out[i] = (E1 - E2)/(aa[i]-bb[i]);
  }
  return out;
}

// [[Rcpp::export]]
mat Hdelta(vec& delta, mat& G, vec& b, mat& X, mat& Z, vec& Y, vec& lY, vec& beta, vec& tau, vec& w, vec& v, int summax, double eps){
  int n = delta.size();
  
  if(n == 1){
    mat H = mat(1, 1);
    vec dpe = delta + eps, dme = delta - eps;
    H(0,0) = (E_llcmp_delta(dpe, G, b, X, Z, Y, lY, beta, tau, w, v, summax) - 
      2.0 * E_llcmp_delta(delta, G, b, X, Z, Y, lY, beta, tau, w, v, summax) + 
      E_llcmp_delta(dme, G, b, X, Z, Y, lY, beta, tau, w, v, summax))/pow(eps, 2.0);
    return H;
  }else{
    vec e = vec(n);
    e.fill(eps);
    mat H = zeros<mat>(n, n);
    mat ee = diagmat(e);
    // Begin loop
    for(int i = 0; i < n; i++){
      vec ei = ee.col(i);
      vec dpe = delta + ei, dme = delta - ei;
      // Diagonal terms
      H(i, i) = (E_llcmp_delta(dme, G, b, X, Z, Y, lY, beta, tau, w, v, summax) - 
        2.0 * E_llcmp_delta(delta, G, b, X, Z, Y, lY, beta, tau, w, v, summax) + 
        E_llcmp_delta(dpe, G, b, X, Z, Y, lY, beta, tau, w, v, summax))/pow(eps, 2.0);
      for(int j = (i + 1); j < n; j++){
        vec ej = ee.col(j); // Off-diagonals
        vec pp = delta + ei + ej, pm = delta + ei - ej, mp = delta - ei + ej, mm = delta - ei - ej;
        H(i,j) =  (E_llcmp_delta(pp, G, b, X, Z, Y, lY, beta, tau, w, v, summax) - 
          E_llcmp_delta(pm, G, b, X, Z, Y, lY, beta, tau, w, v, summax) - 
          E_llcmp_delta(mp, G, b, X, Z, Y, lY, beta, tau, w, v, summax) + 
          E_llcmp_delta(mm, G, b, X, Z, Y, lY, beta, tau, w, v, summax))/(4.0 * pow(eps, 2.0));
        H(j,i) = H(i,j);
      }
    }
    return H;
  }

}

// The hessian at point estimate for b
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

// misc testing
// [[Rcpp::export]]
vec lspacetest(int start, int end, int num){
  return linspace(start, end, num);
}

/* *******
 * END ***
 * *******/

