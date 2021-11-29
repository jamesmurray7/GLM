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
	vec etazi0 = etazi.elem(y0);
	vec eta0 = eta.elem(y0);
	vec out = -exp(etazi) / (1 + exp(etazi));
	out.elem(y0) = out.elem(y0) + (exp(etazi0 + exp(eta0)))/((exp(etazi0 + exp(eta0))) + 1);
	//vec eta0 = eta.elem(y0); 
	//vec eta1 = eta.elem(y1);
	//vec etazi0 = etazi.elem(y0);
	//vec etazi1 = etazi.elem(y1);
	//out.elem(y0) = exp(etazi0) / (exp(etazi0)+exp(-exp(eta0))) - exp(etazi0) / (1 + exp(etazi0)) ;
	//out.elem(y1) = -exp(etazi1) / (1 + exp(etazi1));
	return out;
}


// Functions for items associated with working out things for b
// [[Rcpp::export]]
double b_logdensity(vec& b, vec& Y, mat& X, mat& Z, mat& Xz, mat& Zz,
					vec& beta, vec& alpha, mat& D, int indzi){
	vec eta = X * beta + Z * b(indzi - 2);
	vec etazi = Xz * alpha + Zz * b(indzi - 1);
	double lhs = as_scalar(sum(zip_logdensity(Y, eta, etazi)));
	double rhs = as_scalar(0.5 * b.t() * D.i() * b);
	return -1.0 * (lhs - rhs);
} 

// The same but with c(beta, alpha) as a variable...
// [[Rcpp::export]]
double b_logdensity2(vec& betaalpha, vec& b, vec& Y, mat& X, mat& Z, mat& Xz, mat& Zz,
				     int beta_length, int alpha_length, mat& D, int indzi){
	vec beta = betaalpha.head(beta_length);
	vec alpha = betaalpha.tail(alpha_length);
	vec eta = X * beta + Z * b(indzi - 2);
	vec etazi = Xz * alpha + Zz * b(indzi - 1);
	double lhs = as_scalar(sum(zip_logdensity(Y, eta, etazi)));
	double rhs = as_scalar(0.5 * b.t() * D.i() * b);
	return -1.0 * (lhs - rhs);
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

// Joint likelihood 
// [[Rcpp::export]]
double joint_density(vec& b, vec& Y, mat& X, mat& Z, mat& Xz, mat& Zz,
					 vec& beta, vec& alpha, mat& D, int indzi,
					 double gamma, rowvec& K, vec& eta, rowvec& haz, mat& KK, mat& Fu, double l0i, int Delta){ // survival stuff
	double lhs = b_logdensity(b, Y, X, Z, Xz, Zz, beta, alpha, D, indzi);
	double temp = 0.0;
	if(Delta == 1) temp = log(l0i);
	double surv = as_scalar(temp + Delta * (K * eta + sum(gamma * b)) - haz * (exp(KK * eta + Fu * (gamma * b))));
	return lhs - surv;
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
vec Sgammaeta(vec& gammaeta, vec& b, vec& Y, mat& X, mat& Z, mat& Xz, mat& Zz,               
			  vec& beta, vec& alpha, mat& D, int indzi,                                                  // ZIP 
			  rowvec& K, rowvec& haz, mat& KK, mat& Fu, double l0i, int Delta,   // Survival
			  vec& tau, vec& w, vec& v){                                                                 // Quadrature
   // initialise
   vec out = vec(gammaeta.size());
   int gh = w.size();
   double score_gamma = 0.0;
   vec score_eta = vec(gammaeta.size() - 1.0);
   // Extract gamma and eta
   double gamma = as_scalar(gammaeta.head(1));
   vec eta = gammaeta.tail(2);
   // tau items
   uvec zzz = find(tau == 0.0);
   vec tau2 = pow(pow(gamma, 2.0) * tau, -0.5);
   tau2.elem(zzz).zeros();
   // Loop over quadrature nodes
	for(int l = 0; l < gh; l++){
		vec xi = haz.t() % exp(KK * eta + gamma * Fu * b + pow(pow(gamma, 2.0) * tau, 0.5) * v[l]);
		score_gamma += as_scalar(w[l] * xi.t() * (Fu * b) + gamma * w[l] * v[l] * (xi % tau2).t() * tau);
		score_eta += w[l] * KK.t() * xi;
	}
	out.head(1) = Delta * sum(b) - score_gamma;
	out.tail(2) = Delta * K.t() - score_eta;
	return out;
}
// [[Rcpp::export]]
mat Hgammaeta(vec& gammaeta, vec& b, vec& Y, mat& X, mat& Z, mat& Xz, mat& Zz,               
			  vec& beta, vec& alpha, mat& D, int indzi,   
			  rowvec& K, rowvec& haz, mat& KK, mat& Fu, double l0i, int Delta,   
			  vec& tau, vec& w, vec& v, double eps){
	int n = gammaeta.size();
	mat out = zeros<mat>(n, n);
	vec f0 = Sgammaeta(gammaeta, b, Y, X, Z, Xz, Zz, beta, alpha, D, indzi, K, 
	                   haz, KK, Fu, l0i, Delta, tau, w, v);
	for(int i = 0; i < n; i++){
		vec ge = gammaeta;
		double xi = std::max(gammaeta[i], 1.0);
		ge[i] = gammaeta[i] + (eps*xi);
		vec fdiff = Sgammaeta(ge, b, Y, X, Z, Xz, Zz, beta, alpha, D, indzi, K, 
	                   haz, KK, Fu, l0i, Delta, tau, w, v) - f0;
	    out.col(i) = fdiff/(ge[i]-gammaeta[i]);
	}
	return 0.5 * (out + out.t());
}

// Update for (beta, alpha) ----
// [[Rcpp::export]]
vec Sbetaalpha(vec& betaalpha, vec& b, vec& Y, mat& X, mat& Z, mat& Xz, mat& Zz, 
			   int beta_length, int alpha_length, int& indzi, mat& S, vec& w, vec& v){
	vec out = vec(betaalpha.size());
	int gh = w.size();
	// Make tau terms
	vec tau = sqrt(diagvec(join_rows(Z, Zz) * S * join_rows(Z, Zz).t()));
	// Create dummy stores
	vec beta = betaalpha.head(beta_length);
	vec alpha = betaalpha.tail(alpha_length);
	vec Sbeta = vec(beta.size());
	vec Salpha = vec(alpha.size());
	// loop over quadrature nodes
	for(int l = 0; l < gh; l++){
		vec eta = X * beta + Z * b[indzi-2] + v[l] * tau;
		vec etaz = Xz * alpha + Zz * b[indzi-1] + v[l] * tau;
		Sbeta += w[l] * X.t() * zip_Seta(Y, eta, etaz);
		Salpha += w[l] * Xz.t() * zip_Setazi(Y, eta, etaz);
	}
	out.head(beta.size()) = Sbeta;
	out.tail(alpha.size()) = Salpha;
	return out;
}

// [[Rcpp::export]]
vec S2betaalpha(vec& betaalpha, vec& b, vec& Y, mat& X, mat& Z, mat& Xz, mat& Zz, 
	            int beta_length, int alpha_length, int& indzi, mat& S, vec& w, vec& v){
	vec out = vec(betaalpha.size());
	int gh = w.size();
	uvec y0 = find(Y == 0);
	uvec y1 = find(Y > 0);
	// Extract beta and alpha, and create tau
	vec beta = betaalpha.head(beta_length);
	vec alpha = betaalpha.tail(alpha_length);
	vec Sbeta = vec(beta.size());
	vec Salpha = vec(alpha.size());
	vec tau = sqrt(diagvec(join_rows(Z, Zz) * S * join_rows(Z, Zz).t()));
	for(int l = 0; l < gh; l++){
		vec mu = exp(X * beta + Z * b[indzi - 2] + tau * v[l]);
		vec lambda = exp(Xz * alpha + Zz * b[indzi - 1] + tau * v[l]);
		vec lambda_exp_mu = lambda % trunc_exp(mu);
		lambda_exp_mu.replace(datum::inf, 1e200);
		Sbeta += -w[l] * X.rows(y0).t() * (mu.elem(y0) / (lambda_exp_mu.elem(y0) + 1.0)) + w[l] * X.rows(y1).t() * (Y.elem(y1)-mu.elem(y1));
		Salpha += w[l] * Xz.rows(y0).t() * (lambda_exp_mu.elem(y0) / (lambda_exp_mu.elem(y0) + 1)) - w[l] * Xz.t() * (lambda/(1+lambda));
	}
	out.head(beta.size()) = Sbeta;
	out.tail(alpha.size()) = Salpha;
	return out;
}		

//vec S3betaalpha(vec& betaalpha, vec& b, vec& Y, mat& X, mat& Z, mat& Xz, mat& Zz, 
	            //int beta_length, int alpha_length, int& indzi, mat& S, vec& w, vec& v){
	//vec out = vec(betaalpha.size());
	//int gh = w.size();
	//uvec y0 = find(Y == 0);
	//uvec y1 = find(Y > 0);
	//// Extract beta and alpha, and create tau
	//vec beta = betaalpha.head(beta_length);
	//vec alpha = betaalpha.tail(alpha_length);
	//vec Sbeta = vec(beta.size());
	//vec Sbeta_rhs = vec(beta.size());
	//vec Salpha = vec(alpha.size());
	//vec tau = sqrt(diagvec(join_rows(Z, Zz) * S * join_rows(Z, Zz).t()));
	//for(int l = 0; l < gh; l++){
		//vec eta = X * beta + Z * b[indzi - 2] + tau * v[l];
		//vec etaz = Xz * alpha + Z * b[indzi - 1] + tau * v[l];
		//Sbeta += w[l] * X.rows(y0).t() * (
			//(exp(-exp(eta.elem(y0))) % exp(eta.elem(y0))) / 
			//(exp(-exp(eta.elem(y0))) + exp(etaz.elem(y0)))
		//);
		//Sbeta_rhs += w[l] * X.rows(y1).t() * exp(eta.elem(y1));
		//Salpha += w[l] * Xz.t() * (exp(etaz) / (1 + exp(etaz))) + 
		          //w[l] * Xz.rows(y0).t() * (
					//exp(etaz.elem(y0)) / (exp(etaz.elem(y0)) + exp(-exp(eta.elem(y0))))
		          //);
	//}
	//out.head(beta.size()) = -Sbeta + X.rows(y1).t() * Y.elem(y1) - Sbeta_rhs;
	//out.tail(alpha.size()) = -Salpha;
	//return out;
//}				   
		   

// [[Rcpp::export]]
mat Hbetaalpha(vec& betaalpha, vec& b, vec& Y, mat& X, mat& Z, mat& Xz, mat& Zz, 
			   int beta_length, int alpha_length, int& indzi, mat& S, vec& w, vec& v, double eps){
	int n = betaalpha.size();
	mat out = zeros<mat>(n, n);
	vec f0 = S2betaalpha(betaalpha, b, Y, X, Z, Xz, Zz, beta_length, alpha_length,
						indzi, S, w, v);
	for(int i = 0; i < n; i ++){
		vec ba = betaalpha;
		double xi = std::max(betaalpha[i], 1.0);
		ba[i] = betaalpha[i] + (xi * eps);
		vec fdiff = S2betaalpha(ba, b, Y, X, Z, Xz, Zz, beta_length, alpha_length,
						       indzi, S, w, v) - f0;
		out.col(i) = fdiff / (ba[i] - betaalpha[i]);
	}
	return 0.5 * (out + out.t());
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
		
