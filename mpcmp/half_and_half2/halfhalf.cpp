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
  // vec js = linspace(1, summax, summax);
  double j = 0;
  while(j <= summax){
    double outj = pow(j - mu, 2.0) * pow(lambda, j) / (pow(tgamma(j + 1.0), nu) * Z);
    if(isnan(outj)) break;
    out += outj;
    j++;
  }
  return out;
}

// [[Rcpp::export]]
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
vec lambda_appx(vec& mu, vec& nu, int summax){
  int m = mu.size();
  vec lambda = vec(m);
  vec inits = trunc_exp(nu % log(mu + (nu-1.0)/(2.0 * nu)));
  inits.replace(datum::nan, 1.);
  for(int i = 0; i < m; i++){
    double lower = std::max(inits[i] - 10.0, 1e-6);
    double upper = inits[i] + 10.0;
    double f_lower = mu_lambdaZ_eq(lower, mu[i], nu[i], summax);
    double f_upper = mu_lambdaZ_eq(upper, mu[i], nu[i], summax);
    lambda[i] = zeroin_lambda(lower, upper, f_lower, f_upper, mu[i], nu[i], summax, 1e-4, 100);
  }
  return lambda;
}


// Defining a joint density -----------------------------------------------
static double log2pi = log(2.0 * M_PI); 

// [[Rcpp::export]]
double ll_cmp(vec& loglambda, vec& nu, vec& logZ, vec& Y, vec& lY){
  return sum(Y % loglambda - nu % lY - logZ);
}

// The marginal density (and its gradient) of Y|b,delta...
// [[Rcpp::export]]
double marginal_Y(vec& b, mat& X, vec& Y, vec& lY, mat& Z,     // Data matrices
                  vec& beta, double delta, mat& D, int summax){
  // Define mu, nu
  int mi = Y.size();
  vec mu = exp(X * beta + Z * b);
  vec nu = vec(mi, fill::value(exp(delta)));
  // Calculate lambda and logZ
  vec lambda = lambda_appx(mu, nu, summax);
  vec loglambda = log(lambda);
  vec logZ = logZ_c(loglambda, nu, summax);
  // Calculate loglik CMP
  double ll = ll_cmp(loglambda, nu, logZ, Y, lY);
  int q = b.size();
  return -ll -1.0 * as_scalar(-(double)q/2.0 * log2pi - 0.5 * log(det(D)) - 0.5 * b.t() * D.i() * b);
}

// [[Rcpp::export]]
vec marginal_Y_db(vec& b, mat& X, vec& Y, vec& lY, mat& Z,    // Data matrices
                  vec& beta, double delta, mat& D, int summax){
  // Define mu, nu
  int mi = Y.size();
  vec mu = exp(X * beta + Z * b);
  vec nu = vec(mi, fill::value(exp(delta)));
  
  // Calculate lambda, logZ and V
  vec lambda = lambda_appx(mu, nu, summax);
  vec loglambda = log(lambda);
  vec logZ = logZ_c(loglambda, nu, summax);
  vec V = calc_V_vec(mu, lambda, nu, logZ, summax);
  
  // Score of CMP and other constituent distns.
  mat lhs_mat = diagmat(mu) * Z;          
  return -1. *  lhs_mat.t() * ((Y-mu) / V) -1.0 * (-D.i() * b); 
}

// The full joint density f(b,Y,T,...) and its gradient wrt. b
// [[Rcpp::export]]
double joint_density(vec& b, mat& X, vec& Y, vec& lY, mat& Z,             // Data matrices
                     vec& beta, double delta, mat& D,                       // Longit. + RE parameters
                     rowvec& S, mat& SS, rowvec& Fi, mat& Fu, double l0i, rowvec& haz,
                     int Delta, double gamma, vec& zeta, int summax){
  // Define mu, nu
  int mi = Y.size();
  vec mu = exp(X * beta + Z * b);
  vec nu = vec(mi, fill::value(exp(delta)));
  // Rcout << "b: " << b.t() << std::endl;
  // Rcout << "mu: " << mu.t() << std::endl;
  // Rcout << "nu: " << nu.t() << std::endl;
  // Calculate lambda and logZ
  vec lambda = lambda_appx(mu, nu, summax);
  // Rcout << "lam: " << lambda.t() << std::endl;
  vec loglambda = log(lambda);
  vec logZ = logZ_c(loglambda, nu, summax);
  // Rcout << "logZ" << logZ.t() << std::endl;
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
vec joint_density_ddb(vec& b, mat& X, vec& Y, vec& lY, mat& Z,             // Data matrices
                      vec& beta, double delta, mat& D,                       // Longit. + RE parameters
                      rowvec& S, mat& SS, rowvec& Fi, mat& Fu, double l0i, rowvec& haz,
                      int Delta, double gamma, vec& zeta, int summax){
  // Define mu, nu
  int mi = Y.size();
  vec mu = exp(X * beta + Z * b);
  vec nu = vec(mi, fill::value(exp(delta)));

  // Calculate lambda, logZ and V
  vec lambda = lambda_appx(mu, nu, summax);
  vec loglambda = log(lambda);
  vec logZ = logZ_c(loglambda, nu, summax);
  vec V = calc_V_vec(mu, lambda, nu, logZ, summax);
  
  // Score of CMP and other constituent distns.
  mat lhs_mat = diagmat(mu) * Z;          
  vec Score = lhs_mat.t() * ((Y-mu) / V); 
  // Rcout << "CMP Score: " << Score << std::endl;
  
  return -Score + -1.0 * (-D.i() * b + Delta * (Fi.t() * gamma) - 
                          gamma * (Fu.t() * (haz.t() % exp(SS * zeta + Fu * (gamma * b)))));
  
}  

// And take its second derivative via forward differencing ----
// [[Rcpp::export]]
mat joint_density_sdb(vec& b, mat& X, vec& Y, vec& lY, mat& Z,             // Data matrices
                      vec& beta, double delta, mat& D,                       // Longit. + RE parameters
                      rowvec& S, mat& SS, rowvec& Fi, mat& Fu, double l0i, rowvec& haz,
                      int Delta, double gamma, vec& zeta, int summax, double eps){
  int q = b.size();
  mat out = zeros<mat>(q, q);
  // Rcout << "0------------" << std::endl;
  vec f0 = joint_density_ddb(b, X, Y, lY, Z, beta, delta, D, S, SS, Fi, Fu, l0i, haz, Delta, gamma, zeta, summax);
  for(int i = 0; i < q; i++){
    // Rcout << q << "-------------" << std::endl;
    vec bb = b;
    double xi = std::max(1.0, b[i]);
    bb[i] = b[i] + (xi * eps);
    vec fdiff = joint_density_ddb(bb, X, Y, lY, Z, beta, delta, D, S, SS, Fi, Fu, l0i, haz, Delta, gamma, zeta, summax) - f0;
    out.col(i) = fdiff/(bb[i]-b[i]);
  }
  return 0.5 * (out + out.t()); // Ensure symmetry
}

// Forward differencing on objective fn (joint_density) -------------------
double joint_density2(vec b, mat& X, vec& Y, vec& lY, mat& Z,             // Data matrices
                      vec& beta, double delta, mat& D,                       // Longit. + RE parameters
                      rowvec& S, mat& SS, rowvec& Fi, mat& Fu, double l0i, rowvec& haz,
                      int Delta, double gamma, vec& zeta, int summax){
  // Define mu, nu
  int mi = Y.size();
  vec mu = exp(X * beta + Z * b);
  vec nu = vec(mi, fill::value(exp(delta)));
  // Calculate lambda and logZ
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

// [[Rcpp::export]]
mat H_joint_density(vec b, mat& X, vec& Y, vec& lY, mat& Z,
                    vec& beta, double delta, mat& D, rowvec& S, mat& SS, rowvec& Fi,
                    mat& Fu, double l0i, rowvec& haz, int Delta, double gamma,
                    vec zeta, int summax, double eps){
  int q = b.size();
  if(q == 1){
    double val = (joint_density2(b + eps, X, Y, lY, Z, beta, delta, D, S, SS, 
                                 Fi, Fu, l0i, haz, Delta, gamma, zeta, summax) - 
             2. * joint_density2(b, X, Y, lY, Z, beta, delta, D, S, SS, 
                                 Fi, Fu, l0i, haz, Delta, gamma, zeta, summax) + 
                  joint_density2(b - eps, X, Y, lY, Z, beta, delta, D, S, SS, 
                                 Fi, Fu, l0i, haz, Delta, gamma, zeta, summax))/pow(eps, 2.);
    mat H = mat(1, 1, fill::value(val));
    return H;
  }else{
    vec e = vec(q, fill::value(eps));
    mat H = zeros<mat>(q, q), ee = diagmat(e);
    // Begin loop
    for(int i = 0; i < q; i++){
      vec ei = ee.col(i);
      // Diagonal terms
      H(i, i) = (joint_density2(b - ei, X, Y, lY, Z, beta, delta, D, S, SS, 
                 Fi, Fu, l0i, haz, Delta, gamma, zeta, summax) - 
           2. * joint_density2(b, X, Y, lY, Z, beta, delta, D, S, SS, 
                               Fi, Fu, l0i, haz, Delta, gamma, zeta, summax) + 
                joint_density2(b + ei, X, Y, lY, Z, beta, delta, D, S, SS, 
                               Fi, Fu, l0i, haz, Delta, gamma, zeta, summax))/pow(eps, 2.);
      for(int j = (i + 1); j < q; j++){
        vec ej = ee.col(j); // Off-diagonals
        H(i,j) =  (joint_density2(b + ei + ej, X, Y, lY, Z, beta, delta, D, S, SS, 
                   Fi, Fu, l0i, haz, Delta, gamma, zeta, summax) - 
                   joint_density2(b + ei - ej, X, Y, lY, Z, beta, delta, D, S, SS, 
                                  Fi, Fu, l0i, haz, Delta, gamma, zeta, summax) - 
                   joint_density2(b - ei + ej, X, Y, lY, Z, beta, delta, D, S, SS, 
                                  Fi, Fu, l0i, haz, Delta, gamma, zeta, summax) + 
                   joint_density2(b - ei - ej, X, Y, lY, Z, beta, delta, D, S, SS, 
                                  Fi, Fu, l0i, haz, Delta, gamma, zeta, summax))/(4.0 * pow(eps, 2.));
        H(j,i) = H(i,j);
      }
    }
    return H;
  }
}

// Updates for the poisson process parameters -----------------------------
// Conditional expectation on beta
double E_llcmp_beta(vec beta,
                    vec& b, mat& X, mat& Z, vec& Y, 
                    vec& lY, double delta, vec& tau, vec& w, vec& v, int summax){
  // nu unchanging
  int mi = Y.size();
  vec eta = X * beta + Z * b;
  vec nu = vec(mi, fill::value(exp(delta)));
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
vec Sbeta_cdiff(vec beta, vec& b, mat& X, mat& Z, vec& Y, vec& lY, double delta, vec& tau, vec& w, vec& v, int summax, double eps){
  int a = beta.size();
  vec out = vec(a);
  for(int i = 0; i < a; i++){
    vec aa = beta, bb = beta;
    double xi = std::max(std::fabs(beta[i]), 1.0);
    aa[i] = beta[i] + eps * xi;
    bb[i] = beta[i] - eps * xi;
    double E1 = E_llcmp_beta(aa, b, X, Z, Y, lY, delta, tau, w, v, summax);
    double E2 = E_llcmp_beta(bb, b, X, Z, Y, lY, delta, tau, w, v, summax);
    out[i] = (E1 - E2)/(aa[i]-bb[i]);
  }
  return out;
}

// The Hessian via forward differencing
// [[Rcpp::export]]
mat Hbeta(vec beta, vec& b, mat& X, mat& Z, vec& Y, vec& lY, double delta, vec& tau, vec& w, vec& v, int summax, double eps){
  int n = beta.size();
  
  if(n == 1){
    double val = (E_llcmp_beta(beta + eps, b, X, Z, Y, lY, delta, tau, w, v, summax) - 
      2.0 * E_llcmp_beta(beta, b, X, Z, Y, lY, delta, tau, w, v, summax) + 
      E_llcmp_beta(beta - eps, b, X, Z, Y, lY, delta, tau, w, v, summax))/pow(eps, 2.0);
    mat H = mat(1, 1, fill::value(val));
    return H;
  }else{
    vec e = vec(n, fill::value(eps));
    mat H = zeros<mat>(n, n), ee = diagmat(e);
    // Begin loop
    for(int i = 0; i < n; i++){
      vec ei = ee.col(i);
      // Diagonal terms
      H(i, i) = (E_llcmp_beta(beta - ei, b, X, Z, Y, lY, delta, tau, w, v, summax) - 
        2.0 * E_llcmp_beta(beta, b, X, Z, Y, lY, delta, tau, w, v, summax) + 
        E_llcmp_beta(beta + ei, b, X, Z, Y, lY, delta, tau, w, v, summax))/pow(eps, 2.0);
      for(int j = (i + 1); j < n; j++){
        vec ej = ee.col(j); // Off-diagonals
        H(i,j) =  (E_llcmp_beta(beta + ei + ej, b, X, Z, Y, lY, delta, tau, w, v, summax) - 
          E_llcmp_beta(beta + ei - ej, b, X, Z, Y, lY, delta, tau, w, v, summax) - 
          E_llcmp_beta(beta - ei + ej, b, X, Z, Y, lY, delta, tau, w, v, summax) + 
          E_llcmp_beta(beta - ei - ej, b, X, Z, Y, lY, delta, tau, w, v, summax))/(4.0 * pow(eps, 2.0));
        H(j,i) = H(i,j);
      }
    }
    return H;
  }
  
}

// The "raw" score taken with quadrature
// [[Rcpp::export]]
vec Sbeta2(vec& beta, vec& b, mat& X, mat& Z, vec& Y, vec& lY,
    	   double delta, vec& tau, vec& w, vec& v, int summax){
	int P = beta.size(), gh = w.size(), mi = Y.size();
	vec Score = vec(P);
	vec eta = X * beta + Z * b;
	vec nu = vec(mi, fill::value(exp(delta)));
	// Loop over quadrature nodes.
	for(int l = 0; l < gh; l++){
		vec this_eta = eta + v[l] * tau;
		vec mu = exp(this_eta);
		
		// lambda, logZ, V
		vec lam = lambda_appx(mu, nu, summax);
		vec loglam = log(lam);
		vec logZ = logZ_c(loglam, nu, summax);
		vec V = calc_V_vec(mu, lam, nu, logZ, summax);
		
		// Score at this quadrature node/weight
		mat lhs = diagmat(mu) * X;
		vec rhs = (Y-mu)/V;
		Score += w[l] * (lhs.t() * rhs);
	}
	return Score;
}

// The "raw" second derivative taken with quadrature
// [[Rcpp::export]]
mat Hbeta2(vec& beta, vec& b, mat& X, mat& Z, vec& Y, vec& lY,
  double delta, vec& tau, vec& w, vec& v, int summax){
  int P = beta.size(), gh = w.size(), mi = Y.size();
	mat Hess = zeros<mat>(P, P);
  vec eta = X * beta + Z * b;
  vec nu = vec(mi, fill::value(exp(delta)));
  // Loop over quadrature nodes.
  for(int l = 0; l < gh; l++){
      vec this_eta = eta + v[l] * tau;
      vec mu = exp(this_eta);

      // lambda, logZ, V
      vec lam = lambda_appx(mu, nu, summax);
      vec loglam = log(lam);
      vec logZ = logZ_c(loglam, nu, summax);
      vec V = calc_V_vec(mu, lam, nu, logZ, summax);

      // Score at this quadrature node/weight
      mat lhs = X.t() * diagmat(pow(mu,2.)/V);
   		Hess += w[l] * (-lhs * X);
  }
  return Hess;
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

// Update to the dispersion parameter \delta ------------------------------
// The conditional expectation on \b.
// [[Rcpp::export]]
double E_llcmp_delta(double delta, vec& b, mat& X, mat& Z, vec& Y, vec& lY, vec& beta, vec& tau, vec& w, vec& v, int summax){
  // nu unchanging
  int mi = Y.size();
  vec nu = vec(mi, fill::value(exp(delta)));
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
double Sdelta_cdiff(double delta, vec& b, mat& X, mat& Z, vec& Y, vec& lY, vec& beta, vec& tau, vec& w, vec& v, int summax, double eps){
  double aa = delta + eps, bb = delta - eps;
  double E1 = E_llcmp_delta(aa, b, X, Z, Y, lY, beta, tau, w, v, summax);
  double E2 = E_llcmp_delta(bb, b, X, Z, Y, lY, beta, tau, w, v, summax);
  double out = (E1 - E2)/(aa-bb);
  return out;
}

// [[Rcpp::export]]
double Hdelta(double delta, vec& b, mat& X, mat& Z, vec& Y, vec& lY, vec& beta, vec& tau, vec& w, vec& v, int summax, double eps){
  return  (E_llcmp_delta(delta + eps, b, X, Z, Y, lY, beta, tau, w, v, summax) - 
      2.0 * E_llcmp_delta(delta, b, X, Z, Y, lY, beta, tau, w, v, summax) + 
      E_llcmp_delta(delta - eps, b, X, Z, Y, lY, beta, tau, w, v, summax))/pow(eps, 2.0);
}


// void test_num(IntegerVector inds){
//   for(auto num = inds.begin(); num != inds.end(); num++){
//     Rcout << "num is " << num << std::endl;
//     int current = *num;
//     Rcout << "Corresponds to ind: " << current << std::endl;
//   }
// }

// [[Rcpp::export]]
List delta_update(List delta, List b, List X, List Z, List Y, List lY,
                  vec& beta, List tau, vec& w, vec& v, List summax_list, IntegerVector inds){
  vec new_deltas = vec(delta.size()), scores = vec(delta.size()), hess = vec(delta.size());
  for(auto num = inds.begin(); num != inds.end(); num++){
    int current = *num;
    // Unpack neccesary items from provided lists.
    vec b_i = b[current], tau_i = tau[current], Y_i = Y[current], lY_i = lY[current];
    mat X_i = X[current], Z_i = Z[current];
    int summax = summax_list[current];
    double delta_i = delta[current];
    // Calculate first and second derivs...
    double Sd_i = Sdelta_cdiff(delta_i, b_i, X_i, Z_i, Y_i, lY_i, beta, tau_i, w, v, summax, pow(datum::eps, 1./3.));
    double Hd_i = Hdelta(delta_i, b_i, X_i, Z_i, Y_i, lY_i, beta, tau_i, w, v, summax,pow(datum::eps, 1./4.));
    new_deltas[current] = delta_i - (Sd_i/Hd_i);
    scores[current] = Sd_i;
    hess[current] = Hd_i;
  }
  // return new_deltas;
  return List::create(_["new"] = new_deltas, _["scores"] = scores, _["Hessian"] = hess);
}

// [[Rcpp::export]]
cx_mat mat_sqrt(const mat & M){
  return sqrtmat(M);
}

// Expectation and Variance on observed Y ---------------------------------
vec E_Y(vec& loglam, vec& logZ, vec& nu, int summax){
  vec out = vec(loglam.size());
  for(int j = 1; j <= summax; j++){
    out += exp(log((double)j-1.) + ((double)j - 1.) * loglam - nu * lgamma((double)j) - logZ);
  }
  return out;
}

vec V_Y(vec& loglam, vec& logZ, vec& nu, int summax){
  vec out = vec(loglam.size());
  vec EY = E_Y(loglam, logZ, nu, summax);
  for(int j = 1; j <= summax; j++){
    out += exp(2. * log((double)j-1.) + ((double)j - 1.) * loglam - nu * lgamma((double)j) - logZ);
  }
  return out - pow(EY, 2.);
}

//[[Rcpp::export]]
vec expected_variance(vec& mu, vec& nu, int summax){
  vec lam = lambda_appx(mu, nu, summax);
  vec loglam = log(lam);
  vec logZ = logZ_c(loglam, nu, summax);
  return V_Y(loglam, logZ, nu, summax);
}

/* *******
 * END ***
 * *******/

