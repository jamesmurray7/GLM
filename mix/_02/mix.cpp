#include <RcppArmadillo.h>
#include <math.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// 0. Helper functions ------------------------------------------------
// [[Rcpp::export]]
mat makeV(long double v, int l){
  mat vv = diagmat(v * vec(l, fill::ones));
  return vv;
}

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

// 3. Scores for coefficients \beta ------------------------------------
// All-in-one first attempt, largely taking advantage of same covariate structure for all Y
// and just one vector for beta
// [[Rcpp::export]]
vec Sbeta(vec& beta, mat& X, vec& Y_1, vec& Y_2, vec& Y_3,
          mat& Z, vec& b, vec& vars){
	// _g: Gaussian, _b: Binomial, _c: Count (poisson for now)
	vec beta_g =  beta.subvec(0, 3);
	vec beta_b =  beta.subvec(4, 7);
	vec beta_g2 = beta.subvec(8, 11);
	// Random effects
	vec b_g =  b.subvec(0,1);
	vec b_b =  b.subvec(2,3);
	vec b_g2 = b.subvec(4,5);
	
	// Linear predictors
	vec eta_g =  X * beta_g  + Z * b_g;
	vec eta_b =  X * beta_b  + Z * b_b;
	vec eta_g2 = X * beta_g2 + Z * b_g2;
	
	// Make matrices of variances
	mat V1 = makeV(vars[0], Y_1.size());
	mat V2 = makeV(vars[1], Y_3.size());
	
	// Scores
	vec Score_g =  X.t() * Score_eta_gauss(Y_1, eta_g, V1);
	vec Score_b =  X.t() * Score_eta_binom(Y_2, eta_b);
	vec Score_g2 = X.t() * Score_eta_gauss(Y_3, eta_g2, V2);
	
	vec out = join_cols(Score_g, Score_b, Score_g2);
	return out;
}

// Hessian of \beta, done via forward differencing...
// [[Rcpp::export]]
mat Hbeta(vec& beta, mat& X, vec& Y_1, vec& Y_2, vec& Y_3,
          mat& Z, vec& b, vec& vars, double eps){
	int n = beta.size();
	mat out = zeros<mat>(n, n);
	vec f0 = Sbeta(beta, X, Y_1, Y_2, Y_3, Z, b, vars);
	for(int i = 0; i < n; i++){
		vec bb = beta;
		double xi = std::max(beta[i], 1.0);
		bb[i] = beta[i] + (eps * xi);
		vec fdiff = Sbeta(bb, X, Y_1, Y_2, Y_3, Z, b, vars) - f0;
		out.col(i) = fdiff/(bb[i]-beta[i]);
	}
	return 0.5 * (out + out.t());
}

// 4. Setting out joint log density ------------------------------------
// (Negative ll)
// [[Rcpp::export]]
double joint_density(vec& b, mat& X, mat& Z, vec& beta, vec& vars, mat& D,  // Longit data matrices + (co)Variance
					 vec& Y_1, vec& Y_2, vec& Y_3,                                    // The three longitudinal responses.
					 int Delta, rowvec& K, rowvec& Fi, double l0i,                    // Survival
					 mat& KK, mat& Fu, rowvec& haz, vec& gamma, vec& eta){            // -"-
	int q = b.size();
	// _g: Gaussian, _b: Binomial, _c: Count (poisson for now)
	vec beta_g =  beta.subvec(0, 3);
	vec beta_b =  beta.subvec(4, 7);
	vec beta_g2 = beta.subvec(8, 11);
	vec b_g  = b.subvec(0,1);
	vec b_b  = b.subvec(2,3);
	vec b_g2 = b.subvec(4,5);
	
	// Linear predictors
	vec eta_g  = X * beta_g + Z * b_g  ;
	vec eta_b  = X * beta_b + Z * b_b  ;
	vec eta_g2 = X * beta_g2 + Z * b_g2;
	
	// Make residual variance matrices
	mat V1 = makeV(vars[0], Y_1.size());
	mat V2 = makeV(vars[1], Y_3.size());
	
	// Log densities of responses Y_1->Y_3
	double ll_g = gaussian_ll(Y_1, eta_g, V1);
	double ll_b = binomial_ll(Y_2, eta_b);
	double ll_g2 = gaussian_ll(Y_3, eta_g2, V2); 
	double ll = ll_g + ll_b + ll_g2;
	
	double temp = 0.0;
	if(Delta == 1) temp = log(l0i);
	return -ll + -1.0 * as_scalar(
		-q/2.0 * log(2.0 * M_PI) - 0.5 * log(det(D)) - 0.5 * b.t() * D.i() * b + 
		temp + Delta * (K * eta + Fi * (gamma % b)) - haz * exp(KK * eta + repmat(Fu, 1, 3) * (gamma % b))
	);
}

// First derivative w.r.t b
// [[Rcpp::export]]
vec joint_density_ddb(vec& b, mat& X, mat& Z, vec& beta, vec& vars, mat& D,
					  vec& Y_1, vec& Y_2, vec& Y_3, 
					  int Delta, rowvec& K, rowvec& Fi, double l0i,
				 	  mat& KK, mat& Fu, rowvec& haz, vec& gamma, vec& eta){
	// Combine the three LL scores and then subtract d/db(f(RE)) and d/db)(f(survival))
	// _g: Gaussian, _b: Binomial, _c: Count (poisson for now)
	vec beta_g = beta.subvec(0, 3);
	vec beta_b = beta.subvec(4, 7);
	vec beta_g2 = beta.subvec(8, 11);
	vec b_g = b.subvec(0,1);
	vec b_b = b.subvec(2,3);
	vec b_g2 = b.subvec(4,5);
	
	// Linear predictors
	vec eta_g  = X * beta_g + Z * b_g;
	vec eta_b  = X * beta_b + Z * b_b;
	vec eta_g2 = X * beta_g2 + Z * b_g2;
	
	// Make residual variance matrices
	mat V1 = makeV(vars[0], Y_1.size());
	mat V2 = makeV(vars[1], Y_3.size());
	
	// Scores
	vec Score_g = Z.t() * Score_eta_gauss(Y_1, eta_g, V1);
	vec Score_b = Z.t() * Score_eta_binom(Y_2, eta_b);
	vec Score_g2 = Z.t() * Score_eta_gauss(Y_3, eta_g2, V2);
	
	// score for b per Y
	vec SbY = join_cols(Score_g, Score_b, Score_g2); 
	
	return -SbY + -1.0 * (
		-D.i() * b + Delta * (Fi.t() % gamma) - 
		gamma % (repmat(Fu, 1, 3).t() * (haz.t() % exp(KK * eta + repmat(Fu, 1, 3) * (gamma % b))))
	);
}

// Second derivative w.r.t b, done via forward differencing.
// [[Rcpp::export]]
mat joint_density_sdb(vec& b, mat& X, mat& Z, vec& beta, vec& vars, mat& D,
					  vec& Y_1, vec& Y_2, vec& Y_3, 
					  int Delta, rowvec& K, rowvec& Fi, double l0i,
				 	  mat& KK, mat& Fu, rowvec& haz, vec& gamma, vec& eta,
 					  double eps){
	int n = b.size();
	mat out = zeros<mat>(n, n);
	vec f0 = joint_density_ddb(b, X, Z, beta, vars, D, Y_1, Y_2, Y_3,
							   Delta, K, Fi, l0i, KK, Fu, haz, gamma, eta);
	for(int i = 0; i < n; i++){
		vec bb = b;
		double xi = std::max(b[i], 1.0);
		bb[i] = b[i] + (xi * eps);
		vec fdiff = joint_density_ddb(bb, X, Z, beta, vars, D, Y_1, Y_2, Y_3,
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
