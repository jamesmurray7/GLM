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

// The variance, V --------------------------------------------------------
// [[Rcpp::export]]
double calc_V(double mu, double lambda, double nu, double logZ, int summax){
  double term1 = 0.0, term2 = term1;
  vec js = linspace(0, summax, summax + 1);
  for(int j = 0; j < summax; j++){
    // Mean
    term1 += exp(log(js[j] - 1.0) + (js[j] - 1.0) * log(lambda) - nu * lgamma(js[j]) - logZ);
    // Variance LHS
    term2 += exp(2.0 * log(js[j] - 1.0) + (js[j] - 1.0) * log(lambda) - nu * lgamma(js[j]) - logZ);
  }
  return term2 - pow(term1, 2.0);
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

//[[Rcpp::export]]
mat gen_lambda_mat(int N, int summax){
  vec mus = linspace(1.0/((double)N/10.0), 10.00, N), nus = mus;
  mat out = zeros<mat>(N-1, N);
  for(int m = 0; m < (N-1); m ++){
    for(int n = 0; n < N; n++){
      double init = pow(mus[m] + (nus[n] - 1.0) / (2.0 * nus[n]), nus[n]);
      double lower = std::max(1e-6, init - 10.0);
      double upper = std::max(1e3, init + 10.0);
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
  uvec ge10 = find(mu >= 10.0), lt10 = find(mu < 10);
  int a = mu.size(), b = ge10.size(), c = lt10.size();
  vec out = vec(a);
  
  Rcout << "ge10: " << ge10 << std::endl;
  Rcout << "lt10: " << lt10 << std::endl;
  Rcout << "a: " << a << " b: " << b << " c: " << c << std::endl;
  
  // If there are elements >= 10, fix using the approximaton for lambda...
  if(b > 0){
    vec mu_temp = mu.elem(ge10), nu_temp = nu.elem(ge10);
    out.elem(ge10) = lambda_appx(mu_temp, nu_temp, summax);
    Rcout << out << std::endl;
  }
  
  // Simply take matrix lookups for remaining elements
  if(c > 0){
    vec mu_temp = mu.elem(lt10), nu_temp = nu.elem(lt10);
    vec mu2 = mu_fix(mu_temp, N), nu2 = mu_fix(nu_temp, N);         
    vec m = (mu2/(1./((double)N/10.))) - 1.00, n = (nu2/(1./((double)N/10.))) - 1.00;      
    out.elem(lt10) = mat_lookup(m, n, lambda_mat);
    Rcout << out << std::endl;
  }
  
  return out;
}

// [[Rcpp::export]]
vec get_V(vec& mu, vec& nu, int summax, int N, mat& V_mat){
  return vec(mu.size());
}