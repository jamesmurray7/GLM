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

double logfti(vec& b, vec& delta, rowvec& S, mat& SS, rowvec& Fi, mat& Fu,
              double l0i, rowvec& haz, rowvec& WFi, mat& WFu, 
              int Delta, double gamma, double gamma_disp, vec& zeta){
  double temp = 0.0;
  if(Delta == 1) temp = log(l0i);
  
  return as_scalar(
    temp + Delta * (S * zeta + Fi * (gamma * b) + WFi * (gamma_disp * delta)) - 
      haz * exp(SS * zeta + Fu * (gamma * b) + WFu * (gamma_disp * delta))
  );
}

// The full joint density f(b,Y,T,...) and its gradient wrt. b
// [[Rcpp::export]]
double joint_density(vec& b, mat& X, vec& Y, vec& lY, mat& Z, mat& W,     
                     vec& beta, vec& delta, mat& D,                       
                     rowvec& S, mat& SS, rowvec& Fi, mat& Fu, double l0i, rowvec& haz,
                     rowvec& WFi, mat& WFu,
                     int Delta, double gamma, double gamma_disp, vec& zeta, int summax){
  // Define mu, nu
  vec mu = exp(X * beta + Z * b);
  vec nu = exp(W * delta);

  // Calculate lambda and logZ
  vec lambda = lambda_appx(mu, nu, summax);
  vec loglambda = log(lambda);
  vec logZ = logZ_c(loglambda, nu, summax);
  // Calculate loglik CMP and other constituent distns.
  double ll = ll_cmp(loglambda, nu, logZ, Y, lY);
  double ll_fti = logfti(b, delta, S, SS, Fi, Fu, l0i, haz, WFi, WFu, Delta,
                         gamma, gamma_disp, zeta);
  
  int q = b.size();
  return -ll + -1.0 * as_scalar(-(double)q/2.0 * log2pi - 0.5 * log(det(D)) - 0.5 * b.t() * D.i() * b) + ll_fti;
}

// Its first derivative wrt b ----
// [[Rcpp::export]]
vec joint_density_ddb(vec& b, mat& X, vec& Y, vec& lY, mat& Z, mat& W,
                      vec& beta, vec& delta, mat& D,
                      rowvec& S, mat& SS, rowvec& Fi, mat& Fu, double l0i, rowvec& haz,
                      rowvec& WFi, mat& WFu,
                      int Delta, double gamma, double gamma_disp, vec& zeta, int summax){
  // Define mu, nu
  vec mu = exp(X * beta + Z * b);
  vec nu = exp(W * delta); 

  // Calculate lambda, logZ and V
  vec lambda = lambda_appx(mu, nu, summax);
  vec loglambda = log(lambda);
  vec logZ = logZ_c(loglambda, nu, summax);
  vec V = calc_V_vec(mu, lambda, nu, logZ, summax);
  
  // Score of CMP and other constituent distns.
  mat lhs_mat = diagmat(mu) * Z;          
  vec Score = lhs_mat.t() * ((Y-mu) / V); 
  
  return -Score + -1.0 * (-D.i() * b + Delta * (Fi.t() * gamma) - 
                          gamma * (Fu.t() * (haz.t() % exp(SS * zeta + Fu * (gamma * b) + WFu * (gamma_disp * delta)))));
  
}  

// And take its second derivative via forward differencing ----
// [[Rcpp::export]]
mat joint_density_sdb(vec& b, mat& X, vec& Y, vec& lY, mat& Z, mat& W,     
                      vec& beta, vec& delta, mat& D,                       
                      rowvec& S, mat& SS, rowvec& Fi, mat& Fu, double l0i, rowvec& haz,
                      rowvec& WFi, mat& WFu,
                      int Delta, double gamma, double gamma_disp, vec& zeta, int summax, double eps){
  int q = b.size();
  mat out = zeros<mat>(q, q);
  // Rcout << "0------------" << std::endl;
  vec f0 = joint_density_ddb(b, X, Y, lY, Z, W, beta, delta, D, S, SS, Fi, Fu, l0i, haz, WFi, WFu,
                             Delta, gamma, gamma_disp, zeta, summax);
  for(int i = 0; i < q; i++){
    // Rcout << q << "-------------" << std::endl;
    vec bb = b;
    double xi = std::max(1.0, b[i]);
    bb[i] = b[i] + (xi * eps);
    vec fdiff = joint_density_ddb(bb, X, Y, lY, Z, W, beta, delta, D, S, SS, Fi, Fu, l0i, haz, WFi, WFu,
                                  Delta, gamma, gamma_disp, zeta, summax) - f0;
    out.col(i) = fdiff/(bb[i]-b[i]);
  }
  return 0.5 * (out + out.t()); // Ensure symmetry
}

// Forward differencing on objective fn (joint_density) -------------------
double joint_density2(vec b, mat& X, vec& Y, vec& lY, mat& Z, mat& W,     
                      vec& beta, vec& delta, mat& D,                       
                      rowvec& S, mat& SS, rowvec& Fi, mat& Fu, double l0i, rowvec& haz,
                      rowvec& WFi, mat& WFu,
                      int Delta, double gamma, double gamma_disp, vec& zeta, int summax){
  // Define mu, nu
  vec mu = exp(X * beta + Z * b);
  vec nu = exp(W * delta);
  
  // Calculate lambda and logZ
  vec lambda = lambda_appx(mu, nu, summax);
  vec loglambda = log(lambda);
  vec logZ = logZ_c(loglambda, nu, summax);
  // Calculate loglik CMP and other constituent distns.
  double ll = ll_cmp(loglambda, nu, logZ, Y, lY);
  double ll_fti = logfti(b, delta, S, SS, Fi, Fu, l0i, haz, WFi, WFu, Delta,
                         gamma, gamma_disp, zeta);
  int q = b.size();
  return -ll + -1.0 * as_scalar(-(double)q/2.0 * log2pi - 0.5 * log(det(D)) - 0.5 * b.t() * D.i() * b) + ll_fti;
}

// [[Rcpp::export]]
mat H_joint_density(vec b, mat& X, vec& Y, vec& lY, mat& Z, mat& W,     
                    vec& beta, vec& delta, mat& D,                       
                    rowvec& S, mat& SS, rowvec& Fi, mat& Fu, double l0i, rowvec& haz,
                    rowvec& WFi, mat& WFu,
                    int Delta, double gamma, double gamma_disp, vec& zeta, int summax, double eps){
  int q = b.size();
  if(q == 1){
    double val = (joint_density2(b + eps,X, Y, lY, Z, W, beta, delta, D, S, SS, Fi, Fu, l0i, haz, WFi, WFu,
                                 Delta, gamma, gamma_disp, zeta, summax) - 
             2. * joint_density2(b, X, Y, lY, Z, W, beta, delta, D, S, SS, Fi, Fu, l0i, haz, WFi, WFu,
                                 Delta, gamma, gamma_disp, zeta, summax) + 
                  joint_density2(b - eps, X, Y, lY, Z, W, beta, delta, D, S, SS, Fi, Fu, l0i, haz, WFi, WFu,
                                 Delta, gamma, gamma_disp, zeta, summax))/pow(eps, 2.);
    mat H = mat(1, 1, fill::value(val));
    return H;
  }else{
    vec e = vec(q, fill::value(eps));
    mat H = zeros<mat>(q, q), ee = diagmat(e);
    // Begin loop
    for(int i = 0; i < q; i++){
      vec ei = ee.col(i);
      // Diagonal terms
      H(i, i) = (joint_density2(b - ei, X, Y, lY, Z, W, beta, delta, D, S, SS, Fi, Fu, l0i, haz, WFi, WFu,
                 Delta, gamma, gamma_disp, zeta, summax) - 
           2. * joint_density2(b, X, Y, lY, Z, W, beta, delta, D, S, SS, Fi, Fu, l0i, haz, WFi, WFu,
                               Delta, gamma, gamma_disp, zeta, summax) + 
                joint_density2(b + ei, X, Y, lY, Z, W, beta, delta, D, S, SS, Fi, Fu, l0i, haz, WFi, WFu,
                               Delta, gamma, gamma_disp, zeta, summax))/pow(eps, 2.);
      for(int j = (i + 1); j < q; j++){
        vec ej = ee.col(j); // Off-diagonals
        H(i,j) =  (joint_density2(b + ei + ej, X, Y, lY, Z, W, beta, delta, D, S, SS, Fi, Fu, l0i, haz, WFi, WFu,
                   Delta, gamma, gamma_disp, zeta, summax) - 
                   joint_density2(b + ei - ej, X, Y, lY, Z, W, beta, delta, D, S, SS, Fi, Fu, l0i, haz, WFi, WFu,
                                  Delta, gamma, gamma_disp, zeta, summax) - 
                   joint_density2(b - ei + ej, X, Y, lY, Z, W, beta, delta, D, S, SS, Fi, Fu, l0i, haz, WFi, WFu,
                                  Delta, gamma, gamma_disp, zeta, summax) + 
                   joint_density2(b - ei - ej, X, Y, lY, Z, W, beta, delta, D, S, SS, Fi, Fu, l0i, haz, WFi, WFu,
                                  Delta, gamma, gamma_disp, zeta, summax))/(4.0 * pow(eps, 2.));
        H(j,i) = H(i,j);
      }
    }
    return H;
  }
}

// Updates for the poisson process parameters -----------------------------
// The "raw" score taken with quadrature
// [[Rcpp::export]]
vec Sbeta(vec& beta, vec& b, mat& X, mat& Z, mat& W, vec& Y, vec& lY,
    	   vec& delta, vec& tau, vec& w, vec& v, int summax){
	int P = beta.size(), gh = w.size();
	vec Score = vec(P);
	vec eta = X * beta + Z * b;
	vec nu = exp(W * delta);
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
mat Hbeta(vec& beta, vec& b, mat& X, mat& Z, mat& W, vec& Y, vec& lY,
  vec& delta, vec& tau, vec& w, vec& v, int summax){
  int P = beta.size(), gh = w.size();
	mat Hess = zeros<mat>(P, P);
  vec eta = X * beta + Z * b;
  vec nu = exp(W * delta);
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

// The "raw" gradient taken **without** quadrature
// [[Rcpp::export]]
vec Sbeta_noquad(vec& beta, vec& b, mat& X, mat& Z, mat& W, vec& Y, vec& lY,
	vec& delta, int summax){
	vec eta = X * beta + Z * b, nu = exp(W * delta);
	vec mu = exp(eta);

  // lambda, logZ, V
  vec lam = lambda_appx(mu, nu, summax);
  vec loglam = log(lam);
  vec logZ = logZ_c(loglam, nu, summax);
  vec V = calc_V_vec(mu, lam, nu, logZ, summax);

  mat lhs = diagmat(mu) * X;
  vec rhs = (Y-mu)/V;

  return lhs.t() * rhs;
}

// [[Rcpp::export]]
mat Hbeta_noquad(vec& beta, vec& b, mat& X, mat& Z, mat& W, vec& Y, vec& lY,
                 vec& delta, int summax){
  vec eta = X * beta + Z * b, nu = exp(W * delta);
  vec mu = exp(eta);
  
  // lambda, logZ, V
  vec lam = lambda_appx(mu, nu, summax);
  vec loglam = log(lam);
  vec logZ = logZ_c(loglam, nu, summax);
  vec V = calc_V_vec(mu, lam, nu, logZ, summax);
  
  mat lhs = X.t() * diagmat(square(mu)/V);
  
  return -1. * lhs * X;
}


// Updates for the survival pair (gamma, zeta) ----------------------------
// Define the conditional expectation
double Egammazeta(vec& gammazeta, vec& b, mat& Sigma,
                  rowvec& S, mat& SS, mat& Fu, rowvec& Fi, vec& haz, int Delta, 
                  mat& WFu, rowvec& WFi, vec& delta,
                  vec& w, vec& v){
  double g1 = as_scalar(gammazeta.at(0));              // g1 will always be scalar and the first element
  double g2 = as_scalar(gammazeta.at(1));   
  vec z = gammazeta.subvec(2, gammazeta.size() - 1);
  Rcout << "g1" << g1 << std::endl;
  Rcout << "g2" << g2 << std::endl;
  Rcout << "z"  << z.t() << std::endl;
  // determine tau
  vec tau = pow(g1, 2.0) * diagvec(Fu * Sigma * Fu.t());
  double rhs = 0.0;
  for(int l = 0; l < w.size(); l++){
    rhs += w[l] * as_scalar(
      haz.t() * exp(SS * z + Fu * (g1 * b) + v[l] * pow(tau, 0.5) + WFu * (g2 * delta))
    );
  }
  return as_scalar(Delta * (S * z + Fi * (g1 * b) + WFi * (g2 * delta)) - rhs);
}

// Then take Score AND Hessian via forward differencing.
// [[Rcpp::export]]
vec Sgammazeta(vec& gammazeta, vec& b, mat& Sigma,
               rowvec& S, mat& SS, mat& Fu, rowvec& Fi, vec& haz, int Delta, 
               mat& WFu, rowvec& WFi, vec& delta,
               vec& w, vec& v, long double eps){
  int q = gammazeta.size();
  vec out = vec(q);
  double f0 = Egammazeta(gammazeta, b, Sigma, S, SS, Fu, Fi, haz, Delta, WFu, WFi, delta, w, v);
  for(int i = 0; i < q; i++){
    vec ge = gammazeta;
    double xi = std::max(ge[i], 1.0);
    ge[i] = gammazeta[i] + xi * eps;
    double fdiff = Egammazeta(ge, b, Sigma, S, SS, Fu, Fi, haz, Delta, WFu, WFi, delta, w, v) - f0;
    out[i] = fdiff/(ge[i]-gammazeta[i]);
  }
  return out;
}
// [[Rcpp::export]]
mat Hgammazeta(vec& gammazeta, vec& b, mat& Sigma,
                rowvec& S, mat& SS, mat& Fu, rowvec& Fi, vec& haz, int Delta, 
                mat& WFu, rowvec& WFi, vec& delta,
                vec& w, vec& v, long double eps){
  int q = gammazeta.size();
  mat out = zeros<mat>(q, q);
  vec f0 = Sgammazeta(gammazeta, b, Sigma, S, SS, Fu, Fi, haz, Delta, WFu, WFi, delta, w, v, eps);
  for(int i = 0; i < q; i++){
    vec ge = gammazeta;
    double xi = std::max(ge[i], 1.0);
    ge[i] = gammazeta[i] + xi * eps;
    vec fdiff = Sgammazeta(ge, b, Sigma, S, SS, Fu, Fi, haz, Delta, WFu, WFi, delta, w, v, eps) - f0;
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
// NB: "LHS" of log-likelihood is done in R, below is log-likelihood!
// [[Rcpp::export]]
List delta_update(vec& delta, mat& WFu, vec& WFi, int Delta, 
                  mat& SS, mat& Fu, vec& b, double g1, double g2, vec& zeta, 
                  vec& haz, mat& S, vec& w, vec& v){
  int ww = delta.size(), gh = w.size();
  vec Score_lhs = Delta * g2 * WFi, Score_rhs = vec(ww);
  mat Hessian = zeros<mat>(ww, ww);
  vec tau = pow(g1, 2.) * diagvec(Fu * S * Fu.t());
  vec A = SS * zeta + Fu * (g1 * b) + WFu * (g2 * delta);
  // Loop over quadrature nodes...
  Rcout << "score_lhs" << Score_lhs << std::endl;
  Rcout << "score_rhs" << Score_rhs << std::endl;
  for(int l = 0; l < gh; l++){
    // Score
    Score_rhs += -w[l] * g2 * WFu.t() * (
      haz % exp(A + v[l] * pow(tau, 0.5))
    );
    Rcout << "l: " << l <<  "score_rhs" << Score_rhs.t() << std::endl;
    // Hessian
    Hessian += -w[l] * pow(g2, 2.) * (
      diagmat(haz % exp(A + v[l] * pow(tau, 0.5))) * WFu
    ).t() * WFu;
    Rcout << "l: " << l <<  "Hessian" << Hessian << std::endl;
  }
  // Collect
  vec Score = Score_lhs + Score_rhs;
  
  return List::create(_["Score"] = Score, 
                      _["Hessian"] = Hessian);
}


// Expectation and Variance on observed Y ---------------------------------
vec E_Y(vec& loglam, vec& logZ, vec& nu, int summax){
  vec out = vec(loglam.size());
  for(int j = 1; j <= summax; j++){
    out += exp(log((double)j-1.) + ((double)j - 1.) * loglam - nu * lgamma((double)j) - logZ);
  }
  return out;
}

// [[Rcpp::export]]
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

// [[Rcpp::export]]
double marginal_delta_ll(double delta, vec& mu, vec& Y, vec& lY, int xi){
	int mi = Y.size();
	vec nu = vec(mi, fill::value(exp(delta)));
	vec lam = lambda_appx(mu, nu, xi);
	vec loglam = log(lam);
	vec logZ = logZ_c(loglam, nu, xi);
	return -1. * ll_cmp(loglam, nu, logZ, Y, lY);
}

// [[Rcpp::export]]
double to_minimise(double delta, vec& mu, vec& Y, vec& lY, int xi){
  int mi = Y.size();
  double marg = marginal_delta_ll(delta, mu, Y, lY, xi);
  vec nu = vec(mi, fill::value(exp(delta)));
  vec VY = expected_variance(mu, nu, xi);
  double v = var(Y);
  return marg + (std::fabs(sum(v-VY)));
}

// Expectation of log(Y!)
// [[Rcpp::export]]
vec E_lY(vec& loglam, vec& logZ, vec& nu, int summax){
  int mi = nu.size();
  vec out = vec(mi);
  double j = 1.;
  while(j <= summax){
    out += lgamma(j) * exp((j - 1.) * loglam - nu * lgamma(j) - logZ);
    j++;
  }
  return out;
}

// Expectation of Ylog(Y!)
// [[Rcpp::export]]
vec E_YlY(vec& loglam, vec& logZ, vec& nu, int summax){
  int mi = nu.size();
  vec out = vec(mi);
  double j = 1.;
  while(j <= summax){
    out += exp(log(j - 1.) + log(lgamma(j)) + (j - 1.) * loglam - nu * lgamma(j) - logZ);
    j++;
  }
  return out;
}

// Variance of log(Y!)
// [[Rcpp::export]]
vec V_lY(vec& loglam, vec& logZ, vec& nu, vec& ElY, int summax){
  int mi = nu.size();
  vec out = vec(mi);
  double j = 1.;
  while(j <= summax){
    out += pow(lgamma(j), 2.) * exp((j - 1.) * loglam - nu * lgamma(j) - logZ);
    j++;
  }
  return out - pow(ElY, 2.);
}

// [[Rcpp::export]]
mat getW2(vec& A, vec& VY, vec& VlY, mat& W, vec& nu){
	int mi = VY.size(), w = W.n_cols;
	mat W2 = zeros<mat>(w, w);
	vec nu2 = square(nu), A2 = square(A);
	for(int m = 0; m < mi; m++){
		rowvec Wm = W.row(m);
		W2 += ((-A2.at(m) / VY.at(m) + VlY.at(m)) * nu2.at(m)) * Wm.t() * Wm;
	}
	return W2;
}	
  
/* *******
 * END ***
 * *******/

