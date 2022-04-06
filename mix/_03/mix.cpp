#include <RcppArmadillo.h>
#include <math.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// -----------
// mix/_03/mix.cpp
// Simulation scenario outlined in Rustand et al. (2022; ArXiv 220306256).
// -----------

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
vec Sbeta(vec& beta, mat& X, vec& Y_1, vec& Y_2, vec& Y_3, List Z, List b, mat& V){
	// _g: Gaussian, _b: Binomial, _c: Count (poisson for now)
	vec beta_g = beta.subvec(0, 3);
	vec beta_b = beta.subvec(4, 7);
	vec beta_c = beta.subvec(8, 11);
	vec b_g = b[0];
	vec b_b = b[1];
	vec b_c = b[2];
	mat Z_gc = Z["gc"];
	mat Z_b = Z["b"];
	
	// Linear predictors
	vec eta_g = X * beta_g + Z_gc * b_g;
	vec eta_b = X * beta_b + Z_b * b_b;
	vec eta_c = X * beta_c + Z_gc * b_c;
	
	// Scores
	vec Score_g = X.t() * Score_eta_gauss(Y_1, eta_g, V);
	vec Score_b = X.t() * Score_eta_binom(Y_2, eta_b);
	vec Score_c = X.t() * Score_eta_poiss(Y_3, eta_c);
	
	vec out = join_cols(Score_g, Score_b, Score_c);
	return out;
}

// Hessian of \beta, done via forward differencing...
// [[Rcpp::export]]
mat Hbeta(vec& beta, mat& X, vec& Y_1, vec& Y_2, vec& Y_3,
          List Z, List b, mat& V, double eps){
	int n = beta.size();
	mat out = zeros<mat>(n, n);
	vec f0 = Sbeta(beta, X, Y_1, Y_2, Y_3, Z, b, V);
	for(int i = 0; i < n; i++){
		vec bb = beta;
		double xi = std::max(beta[i], 1.0);
		bb[i] = beta[i] + (eps * xi);
		vec fdiff = Sbeta(bb, X, Y_1, Y_2, Y_3, Z, b, V) - f0;
		out.col(i) = fdiff/(bb[i]-beta[i]);
	}
	return 0.5 * (out + out.t());
}

// With quadrature
// [[Rcpp::export]]
vec Sbeta_quad(vec& beta, mat& X, vec& Y_1, vec& Y_2, vec& Y_3,
               List Z, List b, mat& V, mat& S, vec& w, vec& v, double eps){
  // _g: Gaussian, _b: Binomial, _c: Count (poisson for now)
  vec beta_g = beta.subvec(0, 3);
  vec beta_b = beta.subvec(4, 7);
  vec beta_c = beta.subvec(8, 11);
  vec b_g = b[0];
  vec b_b = b[1];
  vec b_c = b[2];
  mat Z_gc = Z["gc"];
  mat Z_b = Z["b"];
  
  // Linear predictors
  vec eta_g = X * beta_g + Z_gc * b_g;
  vec eta_b = X * beta_b + Z_b * b_b;
  vec eta_c = X * beta_c + Z_gc * b_c;
  
  // Scores
  vec Score_g = X.t() * Score_eta_gauss(Y_1, eta_g, V);
  vec Score_c = X.t() * Score_eta_poiss(Y_3, eta_c);
  
  // Quadrature for binomial sub-model.
  int n_b = beta_b.size();
  vec Score_b = vec(n_b);
  double f0 = binomial_ll_quad(beta_b, X, Y_2, Z_b, b_b, S, w, v); 
  for(int i = 0; i < n_b; i++){
    vec bb = beta_b;
    double xi = std::max(1.0, bb[i]);
    bb[i] = beta_b[i] + xi * eps;
    double fdiff = binomial_ll_quad(bb, X, Y_2, Z_b, b_b, S, w, v) - f0;
    Score_b[i] = fdiff/(bb[i]-beta_b[i]);
  }
  vec out = join_cols(Score_g, Score_b, Score_c);
  return out;
}

// [[Rcpp::export]]
mat Hbeta_quad(vec& beta, mat& X, vec& Y_1, vec& Y_2, vec& Y_3,
               List Z, List b, mat& V, mat S, vec& w, vec& v, double eps){
  int n = beta.size();
  mat out = zeros<mat>(n, n);
  vec f0 = Sbeta_quad(beta, X, Y_1, Y_2, Y_3, Z, b, V, S, w, v, eps);
  for(int i = 0; i < n; i++){
    vec bb = beta;
    double xi = std::max(beta[i], 1.0);
    bb[i] = beta[i] + (eps * xi);
    vec fdiff = Sbeta_quad(bb, X, Y_1, Y_2, Y_3, Z, b, V, S, w, v, eps) - f0;
    out.col(i) = fdiff/(bb[i]-beta[i]);
  }
  return 0.5 * (out + out.t());
}



// 4. Setting out joint log density ------------------------------------
// (Negative ll)
// [[Rcpp::export]]
double joint_density(vec& b, mat& X, List Z, vec& beta, mat& V, mat& D, // Longit data matrices + (co)Variance
					 vec& Y_1, vec& Y_2, vec& Y_3, int Delta, rowvec& Fi, double l0i, mat& Fu, rowvec& haz, vec& gamma){       // Y + Survival
	int q = b.size();
	// _g: Gaussian, _b: Binomial, _c: Count (poisson for now)
	vec beta_g = beta.subvec(0, 3);
	vec beta_b = beta.subvec(4, 7);
	vec beta_c = beta.subvec(8, 11);
	vec b_g = b.subvec(0,1);
	vec b_b = b.subvec(2,2);
	vec b_c = b.subvec(3,4);
	mat Z_gc = Z["gc"];
	mat Z_b = Z["b"];
	// Linear predictors
	vec eta_g = X * beta_g + Z_gc * b_g;
	vec eta_b = X * beta_b + Z_b * b_b;
	vec eta_c = X * beta_c + Z_gc * b_c;
	
	// Log densities of responses Y_1->Y_3
	double ll_g = gaussian_ll(Y_1, eta_g, V);
	double ll_b = binomial_ll(Y_2, eta_b);
	double ll_c = poisson_ll(Y_3, eta_c);
	double ll = ll_g + ll_b + ll_c;
	
	double temp = 0.0;
	if(Delta == 1) temp = log(l0i);
	return -ll + -1.0 * as_scalar(
		-q/2.0 * log(2.0 * M_PI) - 0.5 * log(det(D)) - 0.5 * b.t() * D.i() * b + 
		temp + Delta * (Fi * (gamma % b)) - haz * exp(Fu * (gamma % b))
	);
}

// First derivative w.r.t b
// [[Rcpp::export]]
vec joint_density_ddb(vec& b, mat& X, List Z, vec& beta, mat& V, mat& D,          // Longit data matrices + (co)Variance
					            vec& Y_1, vec& Y_2, vec& Y_3,                               // The three longitudinal responses. 
					            int Delta, rowvec& Fi, double l0i, mat& Fu, rowvec& haz, vec& gamma){ // Survival
	// Combine the three LL scores and then subtract d/db(f(RE)) and d/db)(f(survival)) 
	// _g: Gaussian, _b: Binomial, _c: Count (poisson for now)
	vec beta_g = beta.subvec(0, 3);
	vec beta_b = beta.subvec(4, 7);
	vec beta_c = beta.subvec(8, 11);
	vec b_g = b.subvec(0,1);
	vec b_b = b.subvec(2,2);
	vec b_c = b.subvec(3,4);
	mat Z_gc = Z["gc"];
	mat Z_b = Z["b"];
	// Linear predictors
	vec eta_g = X * beta_g + Z_gc * b_g;
	vec eta_b = X * beta_b + Z_b * b_b;
	vec eta_c = X * beta_c + Z_gc * b_c;
	
	// Scores
	vec Score_g = Z_gc.t() * Score_eta_gauss(Y_1, eta_g, V);
	vec Score_b = Z_b.t() * Score_eta_binom(Y_2, eta_b);
	vec Score_c = Z_gc.t() * Score_eta_poiss(Y_3, eta_c);
	
	vec SbY = join_cols(Score_g, Score_b, Score_c); // score for b per Y
	
	return -SbY + -1.0 * (
		-D.i() * b + Delta * (Fi.t() % gamma) - 
		gamma % (Fu.t() * (haz.t() % exp(Fu * (gamma % b))))
	);
}

// Second derivative w.r.t b, done via forward differencing.
// [[Rcpp::export]]
mat joint_density_sdb(vec& b, mat& X, List Z, vec& beta, mat& V, mat& D,          // Longit data matrices + (co)Variance
					            vec& Y_1, vec& Y_2, vec& Y_3,                               // The three longitudinal responses. --> Matrix and then .col(i) better?
					            int Delta, rowvec& Fi, double l0i,                          // Survival
				 	            mat& Fu, rowvec& haz, vec& gamma, double eps){
	int n = b.size();
	mat out = zeros<mat>(n, n);
	vec f0 = joint_density_ddb(b, X, Z, beta, V, D, Y_1, Y_2, Y_3, Delta, Fi, l0i, Fu, haz, gamma);
	for(int i = 0; i < n; i++){
		vec bb = b;
		double xi = std::max(b[i], 1.0);
		bb[i] = b[i] + (xi * eps);
		vec fdiff = joint_density_ddb(bb, X, Z, beta, V, D, Y_1, Y_2, Y_3, Delta, Fi, l0i, Fu, haz, gamma) - f0;
		out.col(i) = fdiff/(bb[i]-b[i]);
	}
	return 0.5 * (out + out.t());
}

// 5. gamma
// Define the conditional expectation and then take Score and Hessian via forward differencing.
// [[Rcpp::export]]
double Egamma(vec& gamma, mat& bmat, List S, mat& Fu, List Fu_list,
              rowvec& Fi, vec& haz, int Delta, vec& w, vec& v){
  int K = S.size();
  int gh = w.size();
  vec tau = vec(haz.size());
  for(int i = 0; i < K; i++){
    mat Si_k = S[i];
    mat Fu_k = Fu_list[i];
    tau += pow(gamma[i], 2.0) * diagvec(Fu_k * Si_k * Fu_k.t());
  }
  double rhs = 0;
  for(int l = 0; l < gh; l++){
    rhs += w[l] * as_scalar(haz.t() * exp(Fu * (bmat.t() * gamma) + v[l] * pow(tau, 0.5)));
  }
  return as_scalar(Delta * (Fi * (bmat.t() * gamma)) - rhs);
}

// [[Rcpp::export]]
vec Sgamma(vec& gamma, mat& bmat, List S, mat& Fu, List Fu_list,
           rowvec& Fi, vec& haz, int Delta, vec& w, vec& v, long double eps){
  vec out = vec(gamma.size());
  double f0 = Egamma(gamma, bmat, S, Fu, Fu_list, Fi, haz, Delta, w, v);
  for(int i = 0; i < gamma.size(); i++){
    vec g = gamma;
    double xi = std::max(g[i], 1.0);
    g[i] = gamma[i] + (xi * eps);
    double fdiff = Egamma(g, bmat, S, Fu, Fu_list, Fi, haz, Delta, w, v) - f0;
    out[i] = fdiff / (g[i] - gamma[i]);
  }
  return out;
}

// [[Rcpp::export]]
mat Hgamma(vec& gamma, mat& bmat, List S, mat& Fu, List Fu_list,
           rowvec& Fi, vec& haz, int Delta, vec& w, vec& v, long double eps){
  mat out = zeros<mat>(gamma.size(), gamma.size());
  vec f0 = Sgamma(gamma, bmat, S, Fu, Fu_list, Fi, haz, Delta, w, v, eps);
  for(int i = 0; i < gamma.size(); i++){
    vec g = gamma;
    double xi = std::max(g[i], 1.0);
    g[i] = gamma[i] + (xi * eps);
    vec fdiff = Sgamma(g, bmat, S, Fu, Fu_list, Fi, haz, Delta, w, v, eps) - f0;
    out.col(i) = fdiff / (g[i] - gamma[i]);
  }
  return 0.5 * (out + out.t());
}


void test(rowvec& v, uvec& row, uvec& col){
  Rcout << "col: " << v.cols(col) << std::endl;
  Rcout << "row: " << v.rows(row) << std::endl;
  mat test = v.cols(col);
}

// [[Rcpp::export]]
mat lambdaUpdate(List survtime, mat& ft, vec& gamma, List S, List b, List inds, vec& w, vec& v){
  int gh = w.size();
  int n = b.size();
  int K = gamma.size();
  mat store = zeros<mat>(ft.n_rows, n); // Initialise the matrix.
  // Start loop over subjects
  for(int i = 0; i < n; i++){
    vec survtimes_i = survtime[i];      // Survived time indices.
    List Si = S[i];
    List bi = b[i];
    for(int j = 0; j < survtimes_i.size(); j++){
      rowvec Fj = ft.row(j);
      double tau = 0.0, mu = tau;
      for(int k = 0; k < K ; k++){
        uvec ids_k = inds[k];
        mat Sik = Si[k];
        vec bik = bi[k];
        mat Fjk = Fj.cols(ids_k);
        mu += gamma[k] * as_scalar((Fjk * bik));
        tau += pow(gamma[k], 2.0) * as_scalar(Fjk * Sik * Fjk.t());
      }
      for(int l = 0; l < gh; l++){
        store(j, i) += w[l] * as_scalar(exp(mu + pow(tau, 0.5) * v[l]));
      }
    }
  }
  return store;
}

// 6. Update for baseline hazard, lambda -------------------------------
// Function for the update to the baseline hazard, \lambda_0(u)
// // [[Rcpp::export]]
// mat lambdaUpdate(List survtimes, vec& ft, vec& gamma, List S, List b, vec& w, vec& v){
// 	int gh = w.size();
// 	int id = b.size();
// 	mat store = zeros<mat>(ft.size(), id); // Initialise the matrix
//  	// Start loop over i subjects
//  	for(int i = 0; i < id; i++){
//  		vec survtimes_i = survtimes[i];    // This id's survived time indices
//  		List Si = S[i];
//  		List bi = b[i];
//  		for(int j = 0; j < survtimes_i.size(); j++){
//  			rowvec Fst = NumericVector::create(1.0, ft[j]);
//  			double tau = 0.0;
//  			vec rhs = NumericVector::create(0.0, 0.0);
// 			// Loop over the K responses
//  			for(int k = 0; k < bi.size(); k++){
// 				mat Sik = Si[k];
// 				vec bik = bi[k];
// 				tau += as_scalar(pow(gamma[k], 2.0) * Fst * Sik * Fst.t());
// 				rhs += gamma[k] * bik;
// 			}
//  			double mu = as_scalar(exp(Ki * eta + Fst * rhs));
//  			// Rcpp::Rcout << "mu = " << mu << std::endl;
//  			for(int l = 0; l < gh; l++){ // Finally loop over gh nodes
//  				store(j,i) += as_scalar(w[l] * mu * exp(v[l] * sqrt(tau)));
//  			}
//  		}
//  	}
//  	return store;
//  }

// 8 Update for var.e
// // [[Rcpp::export]]
// vec vare_update(vec& Y, mat& X, mat& Z, vec& beta, vec& b, vec& tau, vec& w, vec& v){
//   vec out = vec(2);
//   int gh = w.size();
//   for(int l = 0; l < gh; l++){
//     vec resid = X * beta + Z * b + tau * v[l];
//     out[0] += w[l] * as_scalar((Y-resid).t() * (Y-resid));
//   }
//   out[1] = Y.size();
//   return out;
// }
