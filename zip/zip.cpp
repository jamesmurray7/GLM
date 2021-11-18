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
	out(0) = -1.0 * as_scalar(Z.t() * zip_Seta(Y, eta, etazi));
	out(1) = -1.0 * as_scalar(Zz.t() * zip_Setazi(Y, eta, etazi));
	return out + (b.t() * D.i()).t();
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

// beta

