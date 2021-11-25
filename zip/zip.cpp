#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// Functions for items associated with ZIP density
// [[Rcpp::export]]
vec zip_logdensity(vec& Y, vec& eta, vec& etazi){
	uvec y0 = find(Y == 0);
	uvec y1 = find(Y > 0);
	vec out = vec(Y.size());
	vec eta0 = eta.elem(y0); 
	vec eta1 = eta.elem(y1);
	vec etaz0 = etazi.elem(y0);
	out.elem(y0) = log(exp(etaz0) + exp(-exp(eta0)));
	out.elem(y1) = Y.elem(y1) % eta1 - exp(eta1) - lgamma(Y.elem(y1) + 1);
	return out - log(1 + exp(etazi));
}

// [[Rcpp::export]]
vec zip_Seta(vec& Y, vec& eta, vec& etazi){
	uvec y0 = find(Y == 0);
	uvec y1 = find(Y > 0);
	vec out = vec(Y.size());
	vec eta0 = eta.elem(y0); 
	vec eta1 = eta.elem(y1);
	vec etazi0 = etazi.elem(y0);
	out.elem(y0) = -exp(-exp(eta0)) % exp(eta0) / (exp(-exp(eta0)) + exp(etazi0));
	out.elem(y1) = Y.elem(y1) - exp(eta1);
	return out;
}

// [[Rcpp::export]]
vec zip_Setazi(vec& Y, vec& eta, vec& etazi){     
	uvec y0 = find(Y == 0);
	uvec y1 = find(Y > 0);
	vec out = vec(Y.size());
	vec eta0 = eta.elem(y0); 
	vec eta1 = eta.elem(y1);
	vec etazi0 = etazi.elem(y0);
	vec etazi1 = etazi.elem(y1);
	out.elem(y0) = exp(etazi0) / (exp(etazi0)+exp(-exp(eta0))) - exp(etazi0) / (1 + exp(etazi0)) ;
	out.elem(y1) = -exp(etazi1) / (1 + exp(etazi1));
	return out;
}


// Functions for items associated with working out things for b
// [[Rcpp::export]]
double b_logdensity(vec& b, vec& Y, mat& X, mat& Z, mat& Xz, mat& Zz,
					vec& beta, vec& alpha, mat& D, int indzi){
	vec eta = X * beta + Z * b(indzi - 2);
	vec etazi = Xz * alpha + Zz * b(indzi - 1);
	double lhs = -1.0 * as_scalar(sum(zip_logdensity(Y, eta, etazi)));
	double rhs = as_scalar(0.5 * b.t() * D.i() * b);
	return lhs + rhs;
} 

// [[Rcpp::export]]
vec b_score(vec& b, vec& Y, mat& X, mat& Z, mat& Xz, mat& Zz,
			vec& beta, vec& alpha, mat& D, int indzi){
	vec eta = X * beta + Z * b(indzi - 2);
	vec etazi = Xz * alpha + Zz * b(indzi - 1);
	vec out = vec(b.size());
	out(0) =  as_scalar(Z.t() * zip_Seta(Y, eta, etazi));
	out(1) =  as_scalar(Zz.t() * zip_Setazi(Y, eta, etazi));
	return -1.0 * (out - (b.t() * D.i()).t());
}


// Quadrature for update to \beta and \alpha
// \beta
// [[Rcpp::export]]
vec Sbetal(vec& beta, vec& Y, mat& X, mat& Z, 
		   mat& Xz, mat& Zz, vec& b, mat& D, vec& alpha, int& indzi, List S, double w, double v){
	vec out = vec(beta.size());
	mat S1 = S[0];
	mat S2 = S[1];
	// Make tau terms
	vec tau = sqrt(diagvec(Z * S1 * Z.t()));
	vec tauz = sqrt(diagvec(Zz * S2 * Zz.t()));
	// eta at these quadratures
	vec eta = X * beta + Z * b[indzi-2] + tau * v;
	vec etaz = Xz * alpha + Zz * b[indzi-1] + tauz * v;
	return w * X.t() * zip_Seta(Y, eta, etaz);
}

// \alpha
// [[Rcpp::export]]
vec Salphal(vec& alpha, vec& Y, mat& X, mat& Z, 
		   mat& Xz, mat& Zz, vec& beta, mat& D, vec& b, int& indzi, List S, double w, double v){
	vec out = vec(beta.size());
	mat S1 = S[0];
	mat S2 = S[1];
	// Make tau terms
	vec tau = sqrt(diagvec(Z * S1 * Z.t()));
	vec tauz = sqrt(diagvec(Zz * S2 * Zz.t()));
	// eta at these quadratures
	vec eta = X * beta + Z * b[indzi-2] + tau * v;
	vec etaz = Xz * alpha + Zz * b[indzi-1] + tauz * v;
	return w * Xz.t() * zip_Setazi(Y, eta, etaz);
}

// [[Rcpp::export]]
List Sbeta_alpha(vec& b, vec& Y, mat& X, mat& Z, 
				mat& Xz, mat& Zz, vec& beta, mat& D, vec& alpha, int& indzi, List S, vec& w, vec& v){
	List out;
	mat S1 = S[0];
	mat S2 = S[1];
	int gh = w.size();
	// Make tau terms
	vec tau = sqrt(diagvec(Z * S1 * Z.t()));
	vec tauz = sqrt(diagvec(Zz * S2 * Zz.t()));
	// Create dummy stores
	vec Sbeta = vec(beta.size());
	vec Salpha = vec(alpha.size());
	// loop over quadrature nodes
	for(int l = 0; l < gh; l++){
		vec eta = X * beta + Z * b[indzi-2] + v[l] * tau;
		vec etaz = Xz * alpha + Zz * b[indzi-1] + v[l] * tauz;
		Sbeta += w[l] * X.t() * zip_Seta(Y, eta, etaz);
		Salpha += w[l] * Xz.t() * zip_Setazi(Y, eta, etaz);
	}
	out["Sbeta"] = Sbeta;
	out["Salpha"] = Salpha;
	return out;
}

// Forward differencing for Hessians ---------
// b
// [[Rcpp::export]]
mat fd_b(vec& b, vec& Y, mat& X, mat& Z, mat& Xz, mat& Zz,
		 vec& beta, vec& alpha, mat& D, int indzi, double eps){
	int n = b.size();
	mat out = zeros<mat>(n, n);
	vec f0 = b_score(b, Y, X, Z, Xz, Zz, beta, alpha, D, indzi);
	for(int i = 0; i < n; i++){
	  vec b1 = b;
		double xi = std::max(b[i], 1.0);
		b1[i] = b[i] + (eps*xi);
		vec fdiff = b_score(b1, Y, X, Z, Xz, Zz, beta, alpha, D, indzi) - f0;
		out.col(i) = fdiff/(b1[i]-b[i]);
	}
	return 0.5 * (out + out.t());
} 

// beta @ current quadrature l
// [[Rcpp::export]]
mat fd_beta(vec& beta, vec& b, vec& Y, mat& X, mat& Z, mat& Xz, mat& Zz,
		    vec& alpha, mat& D, List S, int indzi, double w, double v, double eps){
	int n = beta.size();
	mat out = zeros<mat>(n, n);
	vec f0 = Sbetal(beta, Y, X, Z, Xz, Zz, b, D, alpha, indzi, S, w, v);
	for(int i = 0; i < n; i++){
	  vec beta1 = beta;
	  double xi = std::max(beta[i], 1.0);
	  beta1[i] = beta[i] + (eps*xi);
	  vec fdiff = Sbetal(beta1, Y, X, Z, Xz, Zz, b, D, alpha, indzi, S, w, v) - f0;
	  out.col(i) = fdiff/(beta1[i]-beta[i]);
	}
	return 0.5 * (out + out.t());
}

// alpha @ current quadrature l
// [[Rcpp::export]]
mat fd_alpha(vec& alpha, vec& b, vec& Y, mat& X, mat& Z, mat& Xz, mat& Zz,
		     vec& beta, mat& D, List S, int indzi, double w, double v, double eps){
	int n = alpha.size();
	mat out = zeros<mat>(n, n);
	vec f0 = Salphal(alpha, Y, X, Z, Xz, Zz, beta, D, b, indzi, S, w, v);
	for(int i = 0; i < n; i++){
	  vec alpha1 = alpha;
	  double xi = std::max(beta[i], 1.0);
	  alpha1[i] = alpha[i] + (eps*xi);
	  vec fdiff = Salphal(alpha1, Y, X, Z, Xz, Zz, beta, D, b, indzi, S, w, v) - f0;
	  out.col(i) = fdiff/(alpha1[i]-alpha[i]);
	}
	return 0.5 * (out + out.t());
}


// Bring all together (this could be more efficient - tau defined in each separate function but could be done locally)
// [[Rcpp::export]]
List beta_alpha_update(vec& beta, vec& alpha, vec& b, vec& Y, mat& X, mat& Z, mat& Xz,
                  mat& Zz, mat& D, List S, int indzi, vec& w, vec& v, double eps){
	List out;
	int nbeta = beta.size();
	int nalpha = alpha.size();
	int gh = w.size();
	// Create dummy stores
	vec Sbeta = vec(nbeta);
	vec Salpha = vec(nalpha);
	mat Hbeta = zeros<mat>(nbeta, nbeta);
	mat Halpha = zeros<mat>(nalpha, nalpha);
	// Loop over quadrature nodes
	for(int l = 0; l < gh; l++){
		Sbeta += Sbetal(beta, Y, X, Z, Xz, Zz, b, D, alpha, indzi, S, w[l], v[l]);
		Salpha += Salphal(alpha, Y, X, Z, Xz, Zz, beta, D, b, indzi, S, w[l], v[l]);
		Hbeta += fd_beta(beta, b, Y, X, Z, Xz, Zz, alpha, D, S, indzi, w[l], v[l], eps);
		Halpha += fd_alpha(alpha, b, Y, X, Z, Xz, Zz, beta, D, S, indzi, w[l], v[l], eps);
	}
	out["Sbeta"] = Sbeta;
	out["Salpha"] = Salpha;
	out["Hbeta"] = Hbeta;
	out["Halpha"] = Halpha;
	return out;
}

// Joint likelihood 
// [[Rcpp::export]]
double joint_density(vec& b, vec& Y, mat& X, mat& Z, mat& Xz, mat& Zz,
					 vec& beta, vec& alpha, mat& D, int indzi,
					 double gamma, rowvec& K, vec& eta, rowvec& haz, mat& KK, mat& Fu, double l0i, int Delta){ // survival stuff
	double lhs = b_logdensity(b, Y, X, Z, Xz, Zz, beta, alpha, D, indzi);
	double temp = 0.0;
	if(Delta == 1) temp = log(l0i);
	double surv = as_scalar(temp + Delta * (K * eta + sum(gamma * b)) - haz * (exp(KK * eta + Fu * (gamma * b))));
	return lhs + surv;
}

// Joint likelihood - first derivative wrt b
// [[Rcpp::export]]
vec joint_density_ddb(vec& b, vec& Y, mat& X, mat& Z, mat& Xz, mat& Zz,
					 vec& beta, vec& alpha, mat& D, int indzi,
					 double gamma, rowvec& K, vec& eta, rowvec& haz, mat& KK, mat& Fu, double l0i, int Delta){
	vec lhs = b_score(b, Y, X, Z, Xz, Zz, beta, alpha, D, indzi);
	vec ones = vec(b.size(), 1.0);
	// Rcout << "test1" << Delta * gamma * ones << std::endl;
	// Rcout << "test2" << gamma * Fu.t() * (haz.t() % exp(KK * eta + gamma * Fu * b)) << std::endl;
	vec rhs = Delta * gamma * ones - gamma * Fu.t() * (haz.t() % exp(KK * eta + gamma * Fu * b));
	return lhs + -1.0 * rhs;
}

// And the second derivative wrt b (via forward differencing)...
// [[Rcpp::export]]
mat joint_density_sdb(vec& b, vec& Y, mat& X, mat& Z, mat& Xz, mat& Zz,
					 vec& beta, vec& alpha, mat& D, int indzi,
					 double gamma, rowvec& K, vec& eta, rowvec& haz, mat& KK, mat& Fu, double l0i, int Delta, double eps){
	int n = b.size();
	mat out = zeros<mat>(n, n);
	vec f0 = joint_density_ddb(b, Y, X, Z, Xz, Zz, beta, alpha, D, indzi, 
	                           gamma, K, eta, haz, KK, Fu, l0i, Delta);
	for(int i = 0; i < n; i++){
		vec b1 = b;
		double xi = std::max(b[i], 1.0);
		b1[i] = b[i] + (eps * xi);
		vec fdiff = joint_density_ddb(b1, Y, X, Z, Xz, Zz, beta, alpha, D, indzi, 
	                           gamma, K, eta, haz, KK, Fu, l0i, Delta) - f0;
		out.col(i) = fdiff/(b1[i]-b[i]);
	}
	return 0.5 * (out + out.t());
}

// Update for (gamma, eta)
// [[Rcpp::export]]
List gamma_eta_update(vec& b, vec& Y, mat& X, mat& Z, mat& Xz, mat& Zz,               
					 vec& beta, vec& alpha, mat& D, int indzi,                                                  // ZIP 
					 double gamma, rowvec& K, vec& eta, rowvec& haz, mat& KK, mat& Fu, double l0i, int Delta,   // Survival
					 vec& tau, vec& w, vec& v){                                                                 // Quadrature
    List out;
	int gh = w.size();
	// Initialise stores
	double score_gamma = 0.0;
	double H_gamma = 0.0;
	vec score_eta = vec(eta.size());
	mat H_eta = zeros<mat>(eta.size(), eta.size());
	vec cross_term = vec(eta.size());
	// Setting up tau terms
	uvec zzz = find(tau == 0.0);
	vec tau2 = pow(pow(gamma, 2.0) * tau, -0.5);
	vec tau3 = pow(pow(gamma, 2.0) * tau, -1.5);
	tau2.elem(zzz).zeros();
	tau3.elem(zzz).zeros();
	// Loop over quadrature nodes
	for(int l = 0; l < gh; l++){
		vec xi = haz.t() % exp(KK * eta + gamma * Fu * b + pow(pow(gamma, 2.0) * tau, 0.5) * v[l]);
		score_gamma += as_scalar(w[l] * xi.t() * (Fu * b) + gamma * w[l] * v[l] * (xi % tau2).t() * tau);
		score_eta += w[l] * KK.t() * xi;
		H_gamma += as_scalar(
			w[l] * b.t() * Fu.t() * (xi % (Fu * b)) + gamma * w[l] * v[l] * b.t() * Fu.t() * (xi % tau2 % tau) + 
			w[l] * v[l] * (xi % tau2).t() * tau + gamma * w[l] * v[l] * (tau % tau2 % xi).t() * tau + 
			w[l] * v[l] * v[l] * gamma * gamma * (tau % tau2 % xi % tau2).t() * tau - 
			w[l] * v[l] * gamma * gamma * (tau % xi % tau3).t() * tau
		);
		H_eta += w[l] * (diagmat(xi) * KK).t() * KK;
		cross_term += w[l] * KK.t() * (xi % (Fu*b)) + gamma * w[l] * v[l] * KK.t() * (xi % tau % tau2);
	}
	// Populate outputted list
	out["Sgamma"] = Delta * sum(b) - score_gamma;
	out["Seta"] = Delta * K.t() - score_eta;
	out["Hgamma"] = -1.0 * H_gamma;
	out["Heta"] = -1.0 * H_eta;
	out["Hgammaeta"] = -1.0 * cross_term;
	return out;
}
		
// \lambda_0 update for after the E-step (heavily commented as adapting old code)
// [[Rcpp::export]]			 
mat lambdaUpdate(const List survtimes, const vec& ft,
				 const double gamma, const vec& eta, 
				 const List K, List S,
				 const List b,
				 const int id, 
				 const vec& w, const vec& v){
	mat store = zeros<mat>(ft.size(), id); // Initialise the matrix
	int gh = w.size();
	// Start loop over i subjects
	for(int i = 0; i < id; i++){
		vec survtimes_i = survtimes[i];    // This id's survived time indices
		// Rcpp::Rcout << "survtimes_i = " << survtimes_i << std::endl;
		// 'Unpack' the Sigma terms
		mat Si = S[i];   
		// mat S1 = as<mat>(Si[0]);
		// mat S2 = as<mat>(Si[1]);
		vec bi = b[i];
		// vec b1i = bi[0];
		// vec b2i = bi[1];
		// Rcpp::Rcout << "b1i = " << b1i << std::endl; 
		// Rcpp::Rcout << "S1 = " << S1 << std::endl;        
		rowvec Ki = K[i];                   // This id's K
		// Rcpp::Rcout << "K = " << Ki << std::endl; 
		for(int j = 0; j < survtimes_i.size(); j++){
			rowvec Fst  = NumericVector::create(1.0, 1.0);
			double tau = as_scalar(sqrt(
				pow(gamma, 2.0) * Fst * Si * Fst.t() 
			));
			double mu = as_scalar(exp(Ki * eta + gamma * (Fst * bi)));
			// Rcpp::Rcout << "mu = " << mu << std::endl;		  
			for(int k = 0; k < gh; k++){
				store(j,i) += as_scalar(w[k] * mu * exp(v[k] * tau));
			}
		}
	}
	
	return store;
}
		
// For testing against numDeriv::hessian in R
// [[Rcpp::export]]
double joint_density_gamma(double gamma, vec& b, vec& Y, mat& X, mat& Z, mat& Xz, mat& Zz,
					 vec& beta, vec& alpha, mat& D, int indzi,
					  rowvec& K, vec& eta, rowvec& haz, mat& KK, mat& Fu, double l0i, int Delta){ // survival stuff
	double lhs = b_logdensity(b, Y, X, Z, Xz, Zz, beta, alpha, D, indzi);
	double temp = 0.0;
	if(Delta == 1) temp = log(l0i);
	double surv = as_scalar(temp + Delta * (K * eta + sum(gamma * b)) - haz * (exp(KK * eta + Fu * (gamma * b))));
	return lhs + surv;
}

// OLD VERSION
//List gamma_eta_update(vec& b, vec& Y, mat& X, mat& Z, mat& Xz, mat& Zz,               
					 //vec& beta, vec& alpha, mat& D, int indzi,                                                  // ZIP 
					 //double gamma, rowvec& K, vec& eta, rowvec& haz, mat& KK, mat& Fu, double l0i, int Delta,   // Survival
					 //vec& tau, vec& w, vec& v){                                                                 // Quadrature
    //List out;
	//int gh = w.size();
	//// Initialise stores
	//double score_gamma = 0.0;
	//double H_gamma = 0.0;
	//vec score_eta = vec(eta.size());
	//mat H_eta = zeros<mat>(eta.size(), eta.size());
	//vec cross_term = vec(eta.size());
	//// Loop over quadrature nodes
	//for(int l = 0; l < gh; l++){
		//vec xi = haz.t() % exp(KK * eta + gamma * Fu * b + tau * v[l]);
		//score_gamma += as_scalar(w[l] * xi.t() * (Fu * b) + w[l] * v[l] * xi.t() * tau);
		//score_eta += w[l] * KK.t() * xi;
		//H_gamma += as_scalar(
			//w[l] * v[l] * (tau % xi).t() * (Fu * b) + w[l] * v[l] * v[l] * (tau % xi).t() * tau + 
			//w[l] * b.t() * Fu.t() * (xi % (Fu * b)) + w[l] * v[l] * b.t() * Fu.t() * (xi % tau)
		//);
		//H_eta += w[l] * (diagmat(xi) * KK).t() * KK;
		//cross_term += w[l] * v[l] * KK.t() * (tau % xi) + w[l] * KK.t() * ((Fu * b) % xi);
	//}
	//// Populate outputted list
	//out["Sgamma"] = Delta * sum(b) - score_gamma;
	//out["Seta"] = Delta * K.t() - score_eta;
	//out["Hgamma"] = -1.0 * H_gamma;
	//out["Heta"] = -1.0 * H_eta;
	//out["Hgammaeta"] = -1.0 * cross_term;
	//return out;
//}
		
