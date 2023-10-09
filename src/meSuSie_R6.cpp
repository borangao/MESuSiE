#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
using namespace std;
// [[Rcpp::export]]
double loglik_cpp(arma::vec V, const arma::mat& betahat, const arma::mat& shat2,const arma::mat& prior_weight, const int nancestry, arma::uvec diag_index,Rcpp::List config_list){		
	// Calculate X-transpose times y matrix using beta hat and shat2
	mat Xty = betahat/shat2;
	// Initialize a correlation matrix of size `nancestry x nancestry`
	mat cor_mat(nancestry,nancestry); 
	// Extract upper triangle indices of the correlation matrix
	uvec upper_indices = trimatu_ind(size(cor_mat));
	// Set the upper triangle elements of the correlation matrix to vector V
	cor_mat(upper_indices) = V;
	// Convert the correlation matrix to a symmetric matrix
	cor_mat = symmatu(cor_mat);
	// Set the diagonal of the correlation matrix to 1
	cor_mat.diag().ones();
	 // Create a diagonal matrix of standard errors (SE)
	mat se_mat = eye(nancestry,nancestry);
	// Update the diagonal elements of V
	V.elem(diag_index-1) =  sqrt(exp(V.elem(diag_index-1))) ;
	// Update the diagonal of the SE matrix using V
	se_mat.diag() = V.elem(diag_index-1);
	// Calculate V matrix using SE matrix and correlation matrix
	arma::mat V_mat = se_mat*cor_mat*se_mat;
	
	// Initialize a matrix to store log Bayes factors (LBF) for each configuration
	arma::mat lbf_mat(shat2.n_rows,config_list.length());
	// Loop through each configuration and compute the log Bayes factor
 for(int n_config=0;n_config<config_list.length();n_config++){
	 vec config_vec = config_list[n_config];
	 config_vec = config_vec -1;
	 uvec config_vec_index = conv_to< uvec >::from( config_vec );
	 
	 arma::vec lbf_config(shat2.n_rows);
	 lbf_config.zeros(); 
	 vec mu = zeros<vec>(shat2.n_rows);
	 // Case where configuration vector has only one element
	 if(config_vec.n_elem==1){
		 lbf_config = arma::log_normpdf(betahat.col(as_scalar(config_vec)),mu,sqrt(as_scalar(V_mat(config_vec_index,config_vec_index)) + shat2.col(as_scalar(config_vec))))-arma::log_normpdf(betahat.col(as_scalar(config_vec)),mu,sqrt(shat2.col(as_scalar(config_vec))));
	 }else if(config_vec.n_elem!=1){
		 arma::mat V_sub = V_mat(config_vec_index,config_vec_index);
		 arma::mat shat2_sub = shat2.cols(config_vec_index);
		 arma::mat Xty_sub = Xty.cols(config_vec_index);

		 double lbf_part_1,lbf_part_2;
		 mat SIGMA,V_SIGMA,post_var,SIGMA_INV;
		 rowvec Xty_var_sub;
		for(int j=0;j<shat2.n_rows;j++){
			SIGMA_INV = diagmat(shat2_sub.row(j));
			SIGMA = inv(SIGMA_INV);
			V_SIGMA = V_sub*SIGMA;
			post_var = SIGMA_INV-inv(SIGMA+SIGMA*V_SIGMA);
			Xty_var_sub = Xty_sub.row(j)*post_var;	
			lbf_part_1 = 0.5*dot(Xty_var_sub.t(), Xty_sub.row(j));
			lbf_part_2 = 0.5*log(det(V_SIGMA+eye(V_sub.n_cols,V_sub.n_cols)));
			lbf_config(j) = lbf_part_1 - lbf_part_2;
		}	
	 }
	 lbf_mat.col(n_config) = lbf_config;
 }
	// Find the maximum LBF across all configurations
	double maxlbf = lbf_mat.max();
	// Compute weights for each configuration using maxlbf
	arma::mat w_multi = exp(lbf_mat - maxlbf)%prior_weight;
	// Calculate the weighted sum of the weights
	double weighted_sum_w = accu(w_multi);
	// Compute the overall log Bayes factor for the model
	double lbf_model = maxlbf + log(weighted_sum_w);
 
	return(-1.0*lbf_model);
 
}
// [[Rcpp::export]]
SEXP mvlmm_reg(arma::mat betahat,arma::mat shat2, arma::mat V_mat){
	// Calculate X-transpose times y matrix using beta hat and shat2
	mat Xty = betahat/shat2;
	// Obtain the number of ancestries from the shat2 matrix
	double nancestry = shat2.n_cols;
	
	// Initialize vectors and matrices to store results
	vec lbf_multi(shat2.n_rows);
	mat post_mean_wmulti(shat2.n_rows,nancestry);
	mat post_mean2_wmulti(shat2.n_rows,nancestry);

	double lbf_part_1,lbf_part_2;
	mat SIGMA,V_SIGMA,post_var,SIGMA_INV;
	rowvec Xty_var;

	for(int i=0;i<shat2.n_rows;i++){
		SIGMA_INV = diagmat(shat2.row(i));
		SIGMA = inv(SIGMA_INV);
		V_SIGMA = V_mat*SIGMA;
		post_var = SIGMA_INV-inv(SIGMA+SIGMA*V_SIGMA);
		Xty_var = Xty.row(i)*post_var;	
		lbf_part_1 = 0.5*dot(Xty_var.t(),Xty.row(i));
		lbf_part_2 = 0.5*log(det(V_SIGMA+eye(nancestry,nancestry)));
		lbf_multi(i) = lbf_part_1 - lbf_part_2;
		post_mean_wmulti.row(i) = Xty_var;
		post_mean2_wmulti.row(i)=square(Xty_var)+post_var.diag().t();
	}
	List res = List::create(Named("lbf") = lbf_multi , Named("post_mean")=post_mean_wmulti,Named("post_mean2")=post_mean2_wmulti);
  
  return(res);
  
}

