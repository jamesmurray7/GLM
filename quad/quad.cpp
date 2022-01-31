#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// 1. Log-likelihoods for the responses `Y.<k>` ------------------------
// Log-likelihood for Gaussian
// [[Rcpp::export]]
double gaussian_ll(vec& b, mat& X, vec& Y, mat& Z, vec& beta, mat& V){
	int m = Y.size();
	vec eta = X * beta + Z * b;
	return as_scalar(
	-m/2.0 * log(2.0 * M_PI) - 0.5 * log(det(V)) - 0.5 * (Y - eta).t() * V.i() * (Y - eta)
	);
}

// Log-likelihood for REs
// [[Rcpp::export]]
double ranef_ll(vec& b, mat& D){
	int q = b.size();
	return as_scalar(
		-q/2.0 * log(2.0 * M_PI) - 0.5 * log(det(D)) - 0.5 * b.t() * D.i() * b
	);
}

// Log-likelihood for time to event
// [[Rcpp::export]]
double tte_ll(vec& b, rowvec& K, mat& KK, rowvec& Fi, mat& Fu,
		      rowvec& l0u, double l0i, int Delta, vec& gamma, vec& eta){
	double temp = 0.0;
	if(Delta == 1) temp = log(l0i);
	return temp + as_scalar(
		Delta * (K * eta + Fi * (gamma % b)) - 
			l0u * exp(KK * eta + repmat(Fu, 1, 3) * (gamma % b))
	);
}

// 2. Setting out the joint (log) density ------------------------------
// [[Rcpp::export]]
double joint_density(vec& b, mat& X, vec& Y, mat& Z, vec& beta, mat& V,
				     mat& D, rowvec& K, mat& KK, rowvec& Fi, mat& Fu,
		             rowvec& l0u, double l0i, int Delta, vec& gamma, vec& eta){
	double llY = gaussian_ll(b, X, Y, Z, beta, V);   // Ys
	double llR = ranef_ll(b, D);                     // Random Effects
	double llT = tte_ll(b, K, KK, Fi, Fu, l0u, l0i, Delta, gamma, eta); // Time to event
	
	return -1.0 * (llY + llR + llT); // Returning neg log lik. for use with ucminf.
} 

// The gradient w.r.t. b
// [[Rcpp::export]]
vec joint_density_db(vec& b, mat& X, vec& Y, mat& Z, vec& beta, mat& V,
				     mat& D, rowvec& K, mat& KK, rowvec& Fi, mat& Fu,
		             rowvec& l0u, double l0i, int Delta, vec& gamma, vec& eta){
	vec resid = Y - X * beta + Z * b;
	vec Score_Y_db = Z.t() * (V.i() * resid);
	vec Score_R_db = -1.0 * D.i() * b;
	vec Score_T_db = Delta * (Fi.t() * gamma) - gamma % (repmat(Fu, 1, 3).t() * (l0u.t() % exp(KK * eta + repmat(Fu, 1, 3) * (gamma % b))));
	
	return -1.0 * (Score_Y_db + Score_R_db + Score_T_db);
}

// The second derivative (Hessian) w.r.t b, taken via forward differencing.
// [[Rcpp::export]]
mat joint_density_sdb(vec& b, mat& X, vec& Y, mat& Z, vec& beta, mat& V,
				      mat& D, rowvec& K, mat& KK, rowvec& Fi, mat& Fu,
		              rowvec& l0u, double l0i, int Delta, vec& gamma, vec& eta){
	int q = b.size();
	mat out = zeros<mat>(q, q);
	vec f0 = joint_density_db(b, X, Y, Z, beta, V, D, K, KK, Fi, Fu, l0u, l0i,
	                          Delta, gamma, eta);
	for(int i = 0; i < q; i++){
		vec bb = b;
		double xi = std::max(b[i], 1.0);
		bb[i] = b[i] + (xi * 1e-4);
		vec fdiff = joint_density_db(bb, X, Y, Z, beta, V, D, K, KK, Fi, Fu, l0u, l0i,
	                          Delta, gamma, eta) - f0;
	    out.col(i) = fdiff/(bb[i]-b[i]);            
	}
	return 0.5 * (out + out.t());
}

// 3. Setting out components for update to (gamma, eta) ----------------
// Conditional expectation
// [[Rcpp::export]]
double Egammaeta(vec& gammaeta, mat& bmat, List S, 
                  rowvec& K, mat& KK, mat& Fu, rowvec& Fi, vec& haz, int Delta,
                  vec& w, vec& v){
  int nK = S.size();
  int gh = w.size();
  vec g = gammaeta.head(nK);
  vec e = gammaeta.tail(2);
  vec tau = vec(Fu.n_rows);
  for(int i = 0; i < S.size(); i++){
	  mat Si = S[i];
	  tau += pow(g[i], 2.0) * diagvec(Fu * Si * Fu.t());
  }
  double rhs = 0.0;
  for(int l = 0; l < gh; l++){
	  rhs += w[l] * as_scalar(haz.t() * exp(KK * e + Fu * (bmat.t() * g) + v[l] * pow(tau, 0.5)));
  }
  return as_scalar(Delta * (K * e + Fi * (bmat.t() * g)) - rhs);
} 

// Score
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

// Hessian
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
 			rowvec Fst = NumericVector::create(1.0, ft[j], pow(ft[j], 2.0));
 			double tau = 0.0;
 			vec rhs = NumericVector::create(0.0, 0.0, 0.0);
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

