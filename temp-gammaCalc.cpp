#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// log-likelihood on survival part for \gamma 
// [[Rcpp::export]]
double gamma_ll(vec& gamma, const int D, const rowvec& K, const mat& KK,
			    const vec& eta, const vec& gr, const rowvec& rvFi, const vec& b,
			    const rowvec& haz, const mat& Fu, List bsplit,
			    List S, const int nK, const vec& w, const vec& v, const int gh){
	// Define mu and tau, which contain summations over nK responses.
	colvec mu_rhs = vec(Fu.n_rows, fill::zeros);
	mat tau = zeros<mat>(Fu.n_rows, Fu.n_rows);
	for(int k = 0; k < nK; k++){
		vec bk = bsplit[k];
		mat Sk = S[k];
		mu_rhs += gamma[k] * Fu * bk;
		tau += pow(gamma[k], 2.0) * Fu * Sk * Fu.t();
	}
	colvec mu = KK * eta + mu_rhs;
	colvec rhs = vec(Fu.n_rows, fill::zeros);
	for(int l = 0; l < gh; l++){
		rhs += w[l] * exp(mu + v[l] * diagvec(tau));
	}
	return as_scalar(
		D * (K * eta + rvFi * (gr % b)) - haz * rhs
	);
}

// [[Rcpp::export]]
Rcpp::NumericVector subset_range(Rcpp::NumericVector x,
                                 int start = 1, int end = 100) {

  // Use the Range function to create a positional index sequence
  return x[Rcpp::Range(start, end)];
}

// [[Rcpp::export]]
double gammaeta_ll(NumericVector gammaeta, const int D, const rowvec& K, const mat& KK,
			       const vec& gr, const rowvec& rvFi, const vec& b,
			       const rowvec& haz, const mat& Fu, List bsplit,
			       List S, const int nK, const vec& w, const vec& v, const int gh){
	// Initialise gamma and eta vectors
	vec gamma = subset_range(gammaeta, 0, (nK-1));
	vec eta = subset_range(gammaeta, (nK), gammaeta.size()-1);
	
	// Rcout << "gamma: " << gamma <<std::endl;
	// Rcout << "eta: " << eta <<std::endl;
	// Define mu and tau, which contain summations over nK responses.
	colvec mu_rhs = vec(Fu.n_rows, fill::zeros);
	mat tau = zeros<mat>(Fu.n_rows, Fu.n_rows);
	for(int k = 0; k < nK; k++){
		vec bk = bsplit[k];
		mat Sk = S[k];
		mu_rhs += gamma[k] * Fu * bk;
		tau += pow(gamma[k], 2.0) * Fu * Sk * Fu.t();
	}
	colvec mu = KK * eta + mu_rhs;
	colvec rhs = vec(Fu.n_rows, fill::zeros);
	for(int l = 0; l < gh; l++){
		rhs += w[l] * exp(mu + v[l] * diagvec(tau));
	}
	return as_scalar(
		D * (K * eta + rvFi * (gr % b)) - haz * rhs
	);
}

// Split out getxi
// [[Rcpp::export]]
colvec getxi(const rowvec& tausurv, const colvec& musurv, double v, const rowvec& haz){
  rowvec xi = haz % (musurv.t() % exp(v * tausurv));
  return xi.t();
}

// Score for \gamma
// \sum_l w_l*xi'*Fu*b_ik+gamma_k*v_l*w_l*[xi % tau2surv]'*tautilde_ik
// Where tau2surv is (sum_k gamma_k*Fu*Sigma_ik*Fu')^(-0.5)
// [[Rcpp::export]]
vec S2gammacalc(const vec& gamma, 
			    const colvec& musurv, const rowvec& tausurv,                    // tausurv is sqrt(diag(full thing))
			    const rowvec& tau2surv, const mat& tautilde,                    // tau2surv as above, tautilde as diag(Fu*S_ik*Fu')
			    const mat& Fu, const rowvec& Fi, const rowvec& haz, List b,
			    const int D, const vec& w, const vec& v, const int nK, const int gh){
	vec out = vec(nK, fill::zeros);
	for(int k = 0; k < nK; k++){
// 		vec tauk = tautilde.row(k);
		vec bk = b[k];
		double lhs = as_scalar(D * Fi * bk);
		double rhs = 0.0;
		for(int l = 0; l < gh; l++){
			colvec xi = getxi(tausurv, musurv, v[l], haz);
			rhs += as_scalar(
				w[l] * xi.t() * Fu * bk + 
				gamma[k] * w[l] * v[l] * ((xi.t() % tau2surv) * tautilde.row(k).t())
			);
		}
		out[k] = lhs - rhs;
	}
	return out;
}

// [[Rcpp::export]]
vec Sgammacalc(const vec& gamma, const int D, const rowvec& tausurv, const colvec& musurv, 
			   const rowvec& haz, const mat& Fu, const rowvec& Fi, const vec& w, const vec& v,
			   List b, const int nK, const int gh){
	vec out = vec(nK, fill::zeros);
	for(int k = 0; k < nK; k++){
		vec bk = b[k];
		double lhs = as_scalar(D * Fi * bk);
		double rhs = 0.0;
		for(int l = 0; l < gh; l++){
			colvec xi = getxi(tausurv, musurv, v[l], haz);
			rhs += as_scalar(w[l] * xi.t() * Fu * bk + gamma[k] * w[l] * v[l] * (xi.t() % tausurv) * xi);
		}
		out[k] = lhs - rhs;
	}
	return out;		   
}


// [[Rcpp::export]]
mat gamma2Calc(const vec& gamma,const mat& tautilde, const rowvec& tausurv, const rowvec& tau2surv, const colvec& musurv, 
               const vec& w, const vec& v, const mat& Fu, const rowvec& haz, List b, int L, int gh){
	mat M = zeros<mat>(L, L);
	for(int i = 0; i < L; i++){
		rowvec bM = b[i];
		for(int j = 0; j < L; j++){
			rowvec bL = b[j];
			for(int k = 0; k < gh; k++){
				const colvec xi = getxi(tausurv, musurv, v(k), haz);
				// Rcpp::Rcout << "xi: " << xi << std::endl;
				const colvec temp1 = xi % tausurv.t() % xi;
				// Rcpp::Rcout << "xi * tau.surv  * xi: " << temp1 << std::endl;
				const colvec temp2 = temp1 % tau2surv.t();
				// Rcpp::Rcout << "xi * tau.surv * xi * tau2surv: " << temp2 << std::endl;
				const colvec temp3 = xi % xi % tau2surv.t();
				// Rcpp::Rcout << "xi * xi * tau2surv: " << temp3 << std::endl;
				const colvec xitau = xi % tausurv.t();
				
				// Rcpp::Rcout << "Whole first line: " << w(k) * bL* Fu.t() * (xi % (Fu * bM.t())) << std::endl;
				// Rcpp::Rcout << "Whole second line: " << gamma[i] * w[k] * v[k] * bL * Fu.t() * (xi % tau2surv.t() % tautilde.row(i).t()) << std::endl;
				// Rcpp::Rcout << "Whole third line: " << 2 * gamma[j] * w[k] * v[k] * temp1.t() * (Fu * bM.t()) << std::endl;
				// Rcpp::Rcout << "Whole fourth line: " << 2 * gamma[i] * gamma[j] * w[k] * v[k] * v[k] * temp2.t() * tautilde.row(i).t() << std::endl;
				// Rcpp::Rcout << "Whole fifth line: " << gamma[i] * gamma[j] * v[k] * w[k] * temp3.t() * tautilde.row(i).t() << std::endl;
				if(i<j){ // The upper triangle
					M(i,j) += as_scalar(w(k) * bL * Fu.t() * (xi % (Fu * bM.t())) + 
						  gamma[i] * w[k] * v[k] * bL * Fu.t() * (xi % tau2surv.t() % tautilde.row(i).t()) +
						  2.0 * gamma[j] * w[k] * v[k] * temp1.t() * (Fu * bM.t()) + 
						  2.0 * gamma[i] * gamma[j] * w[k] * v[k] * v[k] * temp2.t() * tautilde.row(i).t() + 
						  gamma[i] * gamma[j] * v[k] * w[k] * temp3.t() * tautilde.row(i).t());
				}else if(i == j){ // The diagonal terms
					M(i,j) += as_scalar(w(k) * bL * Fu.t() * (xi % (Fu * bM.t())) + 
						      gamma[i] * w[k] * v[k] * bL * Fu.t() * (xi % tau2surv.t() % tautilde.row(i).t()) + 
						      v[k] * w[k] * xitau.t() * xi + 
						      2.0 * gamma[i] * v[k] * w[k] * temp1.t() * (Fu * bL.t()) + 
						      2.0 * gamma[i] * gamma[i] * w[k] * v[k] * v[k] * temp2.t() * tautilde.row(i).t() + 
						      gamma[i] * gamma[i] * v[k] * w[k] * temp3.t() * tautilde.row(i).t());
				}
			}
		M(j,i) = M(i,j);
		}
	}
	return M;
}

// Seta not really necessary (future work for completeness' sake?)
// Ieta
// [[Rcpp::export]]
mat Ietacalc(const int dim, const rowvec& K, const mat& KK, const rowvec& tausurv,
			 const colvec& musurv, const rowvec& haz, const vec& w, const vec& v, const int gh){
	mat M = zeros<mat>(dim, dim);
	for(int l = 0; l < gh; l++){
		vec xi = getxi(tausurv, musurv, v[l], haz);
		//mat ximat = diagmat(xi); <-- uses more memory => slower
		M += w[l] * (KK.t() * diagmat(xi)) * KK;
	}
	return M;
}

// Second derivates d/dgammadeta
// [[Rcpp::export]]
List Igammaetacalc(const int dim, const mat& KK, const rowvec& tausurv, const mat& tautilde,
				   const colvec& musurv, const rowvec& haz, const mat& Fu, List b, 
				   const vec& gamma, const vec& w, const vec& v, const int nK, const int gh){
	List out(nK);
	for(int k = 0; k < nK; k++){
		vec outk = vec(dim, fill::zeros);
		vec bk = b[k];
		for(int l = 0; l < gh; l++){
			colvec xi = getxi(tausurv, musurv, v[l], haz);
			outk += w[l] * KK.t() * ((Fu * bk) % xi) + 
					 2.0 * w[l] * v[l] * gamma[k] * KK.t() * (xi % tautilde.row(k).t() % xi);
		}
		out[k] = outk;
	}
	return out;
}

// [[Rcpp::export]]
List I2gammaetacalc(const vec& gamma, 
					const colvec& musurv, const rowvec& tausurv,                    // tausurv is sqrt(diag(full thing)); taustar = (sumgkFskF)^(-0.5) (i.e. not squared gammas).
					const rowvec& taustar, const mat& tautilde,                     // tautilde as diag(Fu*S_ik*Fu')
					const mat& Fu, const mat& KK, const rowvec& haz, List b,
					const int dim, const vec& w, const vec& v, const int nK, const int gh){
	List out(nK);
	for(int k = 0; k < nK; k++){
		vec outk = vec(dim, fill::zeros);
		vec bk = b[k];
		for(int l = 0; l < gh; l++){
			colvec xi = getxi(tausurv, musurv, v[l], haz);
			outk += w[l] * KK.t() * ((Fu * bk) % xi) + 
					gamma[k] * w[l] * v[l] * KK.t() * (tautilde.row(k).t() % taustar.t() % xi);
		}
	out[k] = outk;
	}
	return out;
}



