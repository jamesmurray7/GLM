#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// \lambda update, does all AFTER the E-step 
// [[Rcpp::export]]
mat lambdaUpdate(const List survtimes, const vec& ft,
				 const vec& gamma, const vec& eta, 
				 const List K, const List S,
				 const List b,
				 const int id, 
				 const vec& w, const vec& v, const int nodes){
	mat store = zeros<mat>(ft.size(), id); // Initialise the matrix
	// Start loop over i subjects
	for(int i = 0; i < id; i++){
		vec survtimes_i = survtimes[i];    // This id's survived time indices
		// Rcpp::Rcout << "survtimes_i = " << survtimes_i << std::endl;
		// 'Unpack' the Sigma terms
		Rcpp::List Si = S[i];   
		mat S1 = as<mat>(Si[0]);
		mat S2 = as<mat>(Si[1]);
		List bi = b[i];
		vec b1i = bi[0];
		vec b2i = bi[1];
		// Rcpp::Rcout << "b1i = " << b1i << std::endl; 
		// Rcpp::Rcout << "S1 = " << S1 << std::endl;        
		rowvec Ki = K[i];                   // This id's K
		// Rcpp::Rcout << "K = " << Ki << std::endl; 
		for(int j = 0; j < survtimes_i.size(); j++){
			// Rcpp::Rcout << "ft[j] = " << ft[j] << std::endl;
			rowvec Fst2  = NumericVector::create(1.0, ft[j]);
			// Rcpp::Rcout << "Fst = " << Fst << std::endl;
			double tau = as_scalar(sqrt(
				pow(gamma[0], 2.0) * Fst2 * S1 * Fst2.t() + 
				pow(gamma[1], 2.0) * Fst2 * S2 * Fst2.t()
			));
			// Rcpp::Rcout << "tau = " << tau << std::endl;
			// Rcpp::Rcout << "Ki * eta = " << Ki * eta << std::endl;
			// Rcpp::Rcout << "gamma[0] * b1i = " << gamma[0] * b1i << std::endl;
			// Rcpp::Rcout << "RHS = " << gamma[0] * b1i + gamma[1] * b2i + gamma[2] * b3i<< std::endl;
			double mu = as_scalar(exp(Ki * eta + (Fst2 * (gamma[0] * b1i) + 
                         Fst2 * (gamma[1] * b2i)).t()));
			// Rcpp::Rcout << "mu = " << mu << std::endl;								  
			for(int k = 0; k < nodes; k++){
				store(j,i) += as_scalar(w[k] * mu * exp(v[k] * tau));
			}
		}
	}
	
	return store;
}
