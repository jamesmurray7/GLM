#include <RcppArmadillo.h>
#include <math.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// 1. Log-likelihoods for the responses `Y.<k>` ------------------------
// Log-likelihood for Gaussian
// [[Rcpp::export]]
double gaussian_ll(vec& Y, vec& eta, mat& V){
	int m = Y.size();
	return as_scalar(
	-m/2.0 * log(2.0 * M_PI) - 0.5 * log(det(V)) - 0.5 * (Y - eta).t() * V.i() * (Y - eta)
	);
}

// Log-likelihood for binomial
// [[Rcpp::export]]
double binomial_ll(vec& Y, vec& eta){
	vec mu = exp(eta)/(1.0 + exp(eta));
	vec out = vec(mu.size());
	for(int i = 0; i < mu.size(); i++){
		out[i] = R::dbinom(Y[i], 1, mu[i], 1);
	}
	return sum(out);
}

// Log-likelihood for Poisson
// [[Rcpp::export]]
double poisson_ll(vec& Y, vec& eta){
	vec out = vec(Y.size());
	for(int i = 0; i < Y.size(); i++){
		out[i] = R::dpois(Y[i], exp(eta[i]), 1);
	}
	return sum(out);
}

// Log-likelihood for Negative Binomial
// [[Rcpp::export]]
double negbin_ll(vec& Y, vec& eta, double theta){
	vec out = vec(Y.size());
	for(int i = 0; i < Y.size(); i++){
		out[i] = R::dnbinom_mu(Y[i], theta, exp(eta[i]), 1);
	}
	return sum(out);
}

// Conditional expectations of log-likelihoods calculated via quadrature ------
// -- Update 29/03/2022:
//    I can only get this to work for binomial sub-model(s).
// --
vec lfact(vec& v){
  vec out = vec(v.size());
  for(int i = 0; i < v.size(); i++){
    out[i] = lgamma(v[i] + 1.0);
  }
  return out;
}


// Poisson with quadrature
// double poisson_ll_quad(vec& beta, mat& X, vec& Y, mat& Z, vec& b,
//                        mat& S, vec& w, vec& v){
//   vec eta = X * beta + Z * b;
//   vec lfactY = lfact(Y);
//   vec out = vec(Y.size()), rhs = out;
//   // quadrature
//   int gh = w.size();
//   vec tau = sqrt(diagvec(Z * S * Z.t()));
//   for(int l = 0; l < gh; l++){
//     vec this_eta = eta + v[l] * tau;
//     rhs += w[l] * exp(this_eta);
//   }
//   out = Y % eta - rhs - lfactY;
//   return sum(out);
// }

// Binomial with Quadrature
// [[Rcpp::export]]
double binomial_ll_quad(vec& beta, mat& X, vec& Y, mat& Z, vec& b,
                        mat& S, vec& w, vec& v){
  vec eta = X * beta + Z * b;
  vec rhs = vec(eta.size());
  // quadrature
  int gh = w.size();
  vec tau = sqrt(diagvec(Z * S * Z.t()));
  for(int l = 0; l < gh; l ++){
    vec this_eta = eta + v[l] * tau;
    rhs += w[l] * log(1.0 + exp(this_eta));
  }
  vec out = Y % eta - rhs;
  return sum(out);
}

// Negative binomial with quadrature
// double negbin_ll_quad(vec& beta, mat& X, vec& Y, mat& Z, vec& b,
//                      mat& S, vec& w, vec& v, double theta){
//   vec eta = X * beta + Z * b;
//   vec rhs = vec(eta.size());
//   // quadrature
//   int gh = w.size();
//   vec tau = sqrt(diagvec(Z * S * Z.t()));
//   for(int l = 0; l < gh; l ++){
//     vec this_eta = eta + v[l] * tau;
//     rhs += w[l] * log(theta + exp(this_eta));
//   }
//   vec out = Y % eta - (theta + Y) % rhs;
//   return sum(out);
// }

// 2. Scores for linear predictors -------------------------------------
// Gaussian -> binom -> poisson -> negbin
// [[Rcpp::export]]
vec Score_eta_gauss(vec& Y, vec& eta, mat& V){
	return V.i() * (Y - eta);
}
// [[Rcpp::export]]
vec Score_eta_binom(vec& Y, vec& eta){
	return Y - exp(eta) / (exp(eta) + 1.0);
}
// [[Rcpp::export]]
vec Score_eta_poiss(vec& Y, vec& eta){
	return Y - exp(eta);
}
// [[Rcpp::export]]
vec Score_eta_negbin(vec& Y, vec& eta, double theta){
	return (theta * (Y - exp(eta))) / (exp(eta) + theta);
}

// 3. Scores for coefficients \beta ------------------------------------
// All-in-one first attempt, largely taking advantage of same covariate structure for all Y
// and just one vector for beta
// [[Rcpp::export]]
vec Sbeta(vec& beta, mat& X, vec& Y_1, vec& Y_2, vec& Y_3,
          mat& Z, vec& b, mat& V, bool nb, double theta){
	// _g: Gaussian, _b: Binomial, _c: Count (poisson for now)
	vec beta_g = beta.subvec(0, 3);
	vec beta_b = beta.subvec(4, 7);
	vec beta_c = beta.subvec(8, 11);
	vec b_g = b.subvec(0,1);
	vec b_b = b.subvec(2,3);
	vec b_c = b.subvec(4,5);
	
	// Linear predictors
	vec eta_g = X * beta_g + Z * b_g;
	vec eta_b = X * beta_b + Z * b_b;
	vec eta_c = X * beta_c + Z * b_c;
	
	// Scores
	vec Score_g = X.t() * Score_eta_gauss(Y_1, eta_g, V);
	vec Score_b = X.t() * Score_eta_binom(Y_2, eta_b);
	vec Score_c = vec(beta_c.size());
	if(nb == TRUE){
		Score_c += X.t() * Score_eta_negbin(Y_3, eta_c, theta);
	}else{
		Score_c += X.t() * Score_eta_poiss(Y_3, eta_c);
	}
	
	vec out = join_cols(Score_g, Score_b, Score_c);
	return out;
}

// Hessian of \beta, done via forward differencing...
// [[Rcpp::export]]
mat Hbeta(vec& beta, mat& X, vec& Y_1, vec& Y_2, vec& Y_3,
          mat& Z, vec& b, mat& V, bool nb, double theta, double eps){
	int n = beta.size();
	mat out = zeros<mat>(n, n);
	vec f0 = Sbeta(beta, X, Y_1, Y_2, Y_3, Z, b, V, nb, theta);
	for(int i = 0; i < n; i++){
		vec bb = beta;
		double xi = std::max(beta[i], 1.0);
		bb[i] = beta[i] + (eps * xi);
		vec fdiff = Sbeta(bb, X, Y_1, Y_2, Y_3, Z, b, V, nb, theta) - f0;
		out.col(i) = fdiff/(bb[i]-beta[i]);
	}
	return 0.5 * (out + out.t());
}

// With quadrature
// [[Rcpp::export]]
vec Sbeta_quad(vec& beta, mat& X, vec& Y_1, vec& Y_2, vec& Y_3,
               mat& Z, vec& b, mat& V, mat& S, vec& w, vec& v, bool nb, double theta, double eps){
  // _g: Gaussian, _b: Binomial, _c: Count (poisson for now)
  vec beta_g = beta.subvec(0, 3);
  vec beta_b = beta.subvec(4, 7);
  vec beta_c = beta.subvec(8, 11);
  vec b_g = b.subvec(0,1);
  vec b_b = b.subvec(2,3);
  vec b_c = b.subvec(4,5);
  
  // Linear predictors
  vec eta_g = X * beta_g + Z * b_g;
  vec eta_b = X * beta_b + Z * b_b;
  vec eta_c = X * beta_c + Z * b_c;
  
  // Scores
  vec Score_g = X.t() * Score_eta_gauss(Y_1, eta_g, V);
  vec Score_c = vec(beta_c.size());
  if(nb == TRUE){
    Score_c += X.t() * Score_eta_negbin(Y_3, eta_c, theta);
  }else{
    Score_c += X.t() * Score_eta_poiss(Y_3, eta_c);
  }
  
  // Quadrature for binomial sub-model.
  int n_b = beta_b.size();
  vec Score_b = vec(n_b);
  double f0 = binomial_ll_quad(beta_b, X, Y_2, Z, b_b, S, w, v); 
  for(int i = 0; i < n_b; i++){
    vec bb = beta_b;
    double xi = std::max(1.0, bb[i]);
    bb[i] = beta_b[i] + xi * eps;
    double fdiff = binomial_ll_quad(bb, X, Y_2, Z, b_b, S, w, v) - f0;
    Score_b[i] = fdiff/(bb[i]-beta_b[i]);
  }
  vec out = join_cols(Score_g, Score_b, Score_c);
  return out;
}

// [[Rcpp::export]]
mat Hbeta_quad(vec& beta, mat& X, vec& Y_1, vec& Y_2, vec& Y_3,
               mat& Z, vec& b, mat& V, mat S, vec& w, vec& v, bool nb, double theta, double eps){
  int n = beta.size();
  mat out = zeros<mat>(n, n);
  vec f0 = Sbeta_quad(beta, X, Y_1, Y_2, Y_3, Z, b, V, S, w, v, nb, theta, eps);
  for(int i = 0; i < n; i++){
    vec bb = beta;
    double xi = std::max(beta[i], 1.0);
    bb[i] = beta[i] + (eps * xi);
    vec fdiff = Sbeta_quad(bb, X, Y_1, Y_2, Y_3, Z, b, V, S, w, v, nb, theta, eps) - f0;
    out.col(i) = fdiff/(bb[i]-beta[i]);
  }
  return 0.5 * (out + out.t());
}

// 4. Setting out joint log density ------------------------------------
// (Negative ll)
// [[Rcpp::export]]
double joint_density(vec& b, mat& X, mat& Z, vec& beta, mat& V, mat& D,           // Longit data matrices + (co)Variance
					 vec& Y_1, vec& Y_2, vec& Y_3, bool nb, double theta,         // The three longitudinal responses.
					 int Delta, rowvec& K, rowvec& Fi, double l0i,                // Survival
					 mat& KK, mat& Fu, rowvec& haz, vec& gamma, vec& eta){        // -"-
	int q = b.size();
	// _g: Gaussian, _b: Binomial, _c: Count (poisson for now)
	vec beta_g = beta.subvec(0, 3);
	vec beta_b = beta.subvec(4, 7);
	vec beta_c = beta.subvec(8, 11);
	vec b_g = b.subvec(0,1);
	vec b_b = b.subvec(2,3);
	vec b_c = b.subvec(4,5);
	
	// Linear predictors
	vec eta_g = X * beta_g + Z * b_g;
	vec eta_b = X * beta_b + Z * b_b;
	vec eta_c = X * beta_c + Z * b_c;
	
	// Log densities of responses Y_1->Y_3
	double ll_g = gaussian_ll(Y_1, eta_g, V);
	double ll_b = binomial_ll(Y_2, eta_b);
	double ll_c = 0.0;
	if(nb == TRUE){
		ll_c += negbin_ll(Y_3, eta_c, theta);
	}else{
	  ll_c += poisson_ll(Y_3, eta_c);
	}
	double ll = ll_g + ll_b + ll_c;
	
	double temp = 0.0;
	if(Delta == 1) temp = log(l0i);
	return -ll + -1.0 * as_scalar(
		-q/2.0 * log(2.0 * M_PI) - 0.5 * log(det(D)) - 0.5 * b.t() * D.i() * b + 
		temp + Delta * (K * eta + Fi * (gamma % b)) - haz * exp(KK * eta + repmat(Fu,1,3) * (gamma % b))
	);
}

// First derivative w.r.t b
// [[Rcpp::export]]
vec joint_density_ddb(vec& b, mat& X, mat& Z, vec& beta, mat& V, mat& D,          // Longit data matrices + (co)Variance
					  vec& Y_1, vec& Y_2, vec& Y_3, bool nb, double theta,     // The three longitudinal responses. --> Matrix and then .col(i) better?
					  int Delta, rowvec& K, rowvec& Fi, double l0i,               // Survival
				 	  mat& KK, mat& Fu, rowvec& haz, vec& gamma, vec& eta){
	// Combine the three LL scores and then subtract d/db(f(RE)) and d/db)(f(survival))
	// _g: Gaussian, _b: Binomial, _c: Count (poisson for now)
	vec beta_g = beta.subvec(0, 3);
	vec beta_b = beta.subvec(4, 7);
	vec beta_c = beta.subvec(8, 11);
	vec b_g = b.subvec(0,1);
	vec b_b = b.subvec(2,3);
	vec b_c = b.subvec(4,5);
	
	// Linear predictors
	vec eta_g = X * beta_g + Z * b_g;
	vec eta_b = X * beta_b + Z * b_b;
	vec eta_c = X * beta_c + Z * b_c;
	
	// Scores
	vec Score_g = Z.t() * Score_eta_gauss(Y_1, eta_g, V);
	vec Score_b = Z.t() * Score_eta_binom(Y_2, eta_b);
	vec Score_c = vec(b_c.size());
	if(nb == TRUE){
		Score_c += Z.t() * Score_eta_negbin(Y_3, eta_c, theta);
	}else{
		Score_c += Z.t() * Score_eta_poiss(Y_3, eta_c);
	}
	
	vec SbY = join_cols(Score_g, Score_b, Score_c); // score for b per Y
	
	return -SbY + -1.0 * (
		-D.i() * b + Delta * (Fi.t() % gamma) - 
		gamma % (repmat(Fu, 1, 3).t() * (haz.t() % exp(KK * eta + repmat(Fu, 1, 3) * (gamma % b))))
	);
}

// Second derivative w.r.t b, done via forward differencing.
// [[Rcpp::export]]
mat joint_density_sdb(vec& b, mat& X, mat& Z, vec& beta, mat& V, mat& D,          // Longit data matrices + (co)Variance
					  vec& Y_1, vec& Y_2, vec& Y_3, bool nb, double theta,        // The three longitudinal responses. --> Matrix and then .col(i) better?
					  int Delta, rowvec& K, rowvec& Fi, double l0i,               // Survival
				 	  mat& KK, mat& Fu, rowvec& haz, vec& gamma, vec& eta,
 					  double eps){
	int n = b.size();
	mat out = zeros<mat>(n, n);
	vec f0 = joint_density_ddb(b, X, Z, beta, V, D, Y_1, Y_2, Y_3, nb, theta, 
							   Delta, K, Fi, l0i, KK, Fu, haz, gamma, eta);
	for(int i = 0; i < n; i++){
		vec bb = b;
		double xi = std::max(b[i], 1.0);
		bb[i] = b[i] + (xi * eps);
		vec fdiff = joint_density_ddb(bb, X, Z, beta, V, D, Y_1, Y_2, Y_3, nb, theta,
							   Delta, K, Fi, l0i, KK, Fu, haz, gamma, eta) - f0;
		out.col(i) = fdiff/(bb[i]-b[i]);
	}
	return 0.5 * (out + out.t());
}

// 5. (gamma, eta)
// Define the conditional expectation and then take Score AND Hessian via forward differencing
// [[Rcpp::export]]
double Egammaeta(vec& gammaeta, mat& bmat, List S,
				rowvec& K, mat& KK, mat& Fu, rowvec& Fi, vec& haz, int Delta, vec& w, vec& v){
	vec g = gammaeta.subvec(0,2);
	vec e = gammaeta.subvec(3,4);
	vec tau = vec(Fu.n_rows);
	for(int i = 0; i < S.size(); i++){
		mat Si = S[i];
		tau += pow(g[i], 2.0) * diagvec(Fu * Si * Fu.t());
	}
	double rhs = 0.0;
	for(int l = 0; l < w.size(); l++){
		rhs += w[l] * as_scalar(haz.t() * exp(KK * e + Fu * (bmat.t() * g) + v[l] * pow(tau, 0.5)));
	}
	return as_scalar(Delta * (K * e + Fi * (bmat.t() * g)) - rhs);
}
// [[Rcpp::export]]
vec Sgammaeta(vec& gammaeta, mat& bmat, List S,
				rowvec& K, mat& KK, mat& Fu, rowvec& Fi, vec& haz, int Delta, vec& w, vec& v, double eps){
	vec out = vec(gammaeta.size());
	double f0 = Egammaeta(gammaeta, bmat, S, K, KK, Fu, Fi, haz, Delta, w, v);
	for(int i = 0; i < gammaeta.size(); i++){
		vec ge = gammaeta;
		double xi = std::max(ge[i], 1.0);
		ge[i] = gammaeta[i] + xi * eps;
		double fdiff = Egammaeta(ge, bmat, S, K, KK, Fu, Fi, haz, Delta, w, v)-f0;
		out[i] = fdiff/(ge[i]-gammaeta[i]);
	}
	return out;
}
// [[Rcpp::export]]
mat Hgammaeta(vec& gammaeta, mat& bmat, List S,
			  rowvec& K, mat& KK, mat& Fu, rowvec& Fi, vec& haz, int Delta, vec& w, vec& v, double eps){
	mat out = zeros<mat>(gammaeta.size(), gammaeta.size());
	vec f0 = Sgammaeta(gammaeta, bmat, S, K, KK, Fu, Fi, haz, Delta, w, v, eps);
	for(int i = 0; i < gammaeta.size(); i++){
		vec ge = gammaeta;
		double xi = std::max(ge[i], 1.0);
		ge[i] = gammaeta[i] + xi * eps;
		vec fdiff = Sgammaeta(ge, bmat, S, K, KK, Fu, Fi, haz, Delta, w, v, eps) - f0;
		out.col(i) = fdiff/(ge[i]-gammaeta[i]);
	}
	return 0.5 * (out + out.t());
}

// 6. Update for baseline hazard, lambda -------------------------------

// Function for the update to the baseline hazard, \lambda_0(u)
//     (adapted from my implementation in Bernhardt work)...
// [[Rcpp::export]]
 mat lambdaUpdate(List survtimes, vec& ft,
 				  vec& gamma, vec& eta, List K, List S, List b, 
 				  vec& w, vec& v){
	int gh = w.size();
	int id = b.size();
	mat store = zeros<mat>(ft.size(), id); // Initialise the matrix
 	// Start loop over i subjects
 	for(int i = 0; i < id; i++){
 		vec survtimes_i = survtimes[i];    // This id's survived time indices   
 		List Si = S[i];
 		List bi = b[i];
 		rowvec Ki = K[i];                  // Start loop over subject i's j survived times     
 		for(int j = 0; j < survtimes_i.size(); j++){
 			rowvec Fst = NumericVector::create(1.0, ft[j]);
 			double tau = 0.0;
 			vec rhs = NumericVector::create(0.0, 0.0);
			// Loop over the K responses 
 			for(int k = 0; k < bi.size(); k++){
				mat Sik = Si[k];
				vec bik = bi[k];
				tau += as_scalar(pow(gamma[k], 2.0) * Fst * Sik * Fst.t());
				rhs += gamma[k] * bik;
			}
 			double mu = as_scalar(exp(Ki * eta + Fst * rhs));
 			// Rcpp::Rcout << "mu = " << mu << std::endl;								  
 			for(int l = 0; l < gh; l++){ // Finally loop over gh nodes
 				store(j,i) += as_scalar(w[l] * mu * exp(v[l] * sqrt(tau)));
 			}
 		}
 	}
 	return store;
 }

// 7. Update for \theta, in case of Y|. ~ NB(mu, theta) parameterisation.
// [[Rcpp::export]]
double ll_theta_quad(double theta, vec& beta, mat& X, vec& Y, mat& Z, vec& b, mat& S,
				     vec& w, vec& v){
	vec rhs = vec(Y.size());
	vec tau = sqrt(diagvec(Z * S * Z.t()));
	vec eta = X * beta + Z * b;
	for(int l = 0; l < w.size(); l++){
		vec this_eta = eta + v[l] * tau;
		rhs += w[l] * log(exp(this_eta) + theta);
	}
	vec out = vec(Y.size());
	for(int i = 0; i < Y.size(); i++){
	   out[i] = lgamma(theta + Y[i]) - lgamma(theta) + theta * log(theta) - (theta + Y[i]) * rhs[i];
	}
	return sum(out);
}

// [[Rcpp::export]]
double Stheta(double theta, vec& beta, mat& X, vec& Y, mat& Z, vec& b, mat& S,
               vec& w, vec& v, double eps){
	double theta1 = theta + eps;
	double theta2 = theta - eps;
	double l1 = ll_theta_quad(theta1, beta, X, Y, Z, b, S, w, v);
	double l2 = ll_theta_quad(theta2, beta, X, Y, Z, b, S, w, v);
	return (l1-l2)/(theta1-theta2);
}

// [[Rcpp::export]]
double Htheta(double theta, vec& beta, mat& X, vec& Y, mat& Z, vec& b, mat& S,
               vec& w, vec& v, double eps){
	double f0 = Stheta(theta, beta, X, Y, Z, b, S, w, v, eps);
	double xi = eps * (abs(theta) + eps);
	double tt = theta + xi;
	double fdiff = Stheta(tt, beta, X, Y, Z, b, S, w, v, eps) - f0;
	return fdiff/(tt - theta);
}

// 8 Update for var.e
// [[Rcpp::export]]
vec vare_update(vec& Y, mat& X, mat& Z, vec& beta, vec& b, vec& tau, vec& w, vec& v){
  vec out = vec(2);
  int gh = w.size();
  for(int l = 0; l < gh; l++){
    vec resid = X * beta + Z * b + tau * v[l];
    out[0] += w[l] * as_scalar((Y-resid).t() * (Y-resid));
  }
  out[1] = Y.size();
  return out;
}
