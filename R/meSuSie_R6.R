#' @title Multiple Ancestry Sum of Single Effects (MESuSiE) Regression
#' @description Performs a multiple ancestry sum of the single effects
#'   regression of multiple ancestry phenotype Y on X.
#'   That is, this function fits the regression model \deqn{Y = \sum_l X
#'   b_l + e,} where the elements of \eqn{e} are \emph{i.i.d.} normal
#'   with zero mean and variance \code{residual_variance}, and the sum
#'   \eqn{\sum_l b_l} is a vector of L effects to be estimated. We
#'   assumes that each \eqn{b_l} has exactly one non-zero
#'   element. We extend the MESuSiE for summary statistics, take only summary statistics 
#'   and LD correlation matrix as input.
#' @param  R_mat_list A list of length N ancestry with each element being correlation matrix with dimension p*p, the column name of the correlation matrix should match to the order of SNP name in summary_stat_list
#' 
#' @param  summary_stat_list A list of length N ancestry with each element being summary statistics. The minimum requirement of summary statistics contains columns of SNP, Beta, Se, Z, and N. 
#' The order of the SNP should match the order of the correlation matrix. MESuSiE required either marginal z-scores and number of individuals from each ancestry, or marginal effect size (Beta) and standard error (Se) of each SNP to 
#' be provides to reconstruct the sufficient statistics.  
#' 
#' @param  L Maximum number of non-zero effects assumed within the region.
#' 
#' @param prior_variance Can be either (1) a vector of length L, or a
#'   scalar, for scaled prior variance when Y is univariate (which
#'   should then be equivalent to \code{\link[susieR]{susie}}); or (2) a
#'   matrix for a simple multivariate regression; or (3) a mixture prior 
#'   from \code{\link{create_mixture_prior}}.
#' 
#' @param residual_variance The residual variance (defaults to be a size N ancestry vector of 1).
#' 
#' @param prior_weights A vector of length p giving the prior
#'   probability that each element is non-zero.
#'  The default setting is that the prior weights are
#'   the same for all variables.
#'   
#' @param ancestry_weight A vector of length 2^{N ancestry}-1, with each element
#' being the prior values for each causal configuration. Suppose there are two ancestries,
#' the vector will be length 3, with element being p1, p2, and p12 which are the prior probability
#' of being uniquely causal in the first ancestry, second ancestry, and causal in both ancestries. The
#' summation of p1, p2, p12 is one. The default prior for the two ancestry setting is 3/7,3/7, and 1/7 to 
#' encourage the finding of ancestry-specific causal variant
#' 
#' @param estimate_residual_variance When
#'   \code{estimate_residual_variance = TRUE}, the residual variance is
#'   estimated; otherwise it is fixed.
#' 
#' @param estimate_prior_variance When \code{estimate_prior_variance =
#'   TRUE}, the prior variance is estimated; otherwise it is
#'   fixed. Currently \code{estimate_prior_variance = TRUE} only works
#'   for univariate Y, or for multivariate Y when the prior variance is
#'   a matrix).
#'   
#' @param cor_method The method to determine the 95% credible set using correlations among 
#' SNPs in the region. The common used can be minimum/mean/median of absolute value of correlation
#' The default is \code{cor_method = min_abs_corr}, which is the
#' minimum of absolute value of correlation allowed in a credible set across ancestries, a commonly
#' used threshold for genotype data in genetics studies. 
#' 
#' @param cor_threshold The threshold of min_abs_corr, defualt is 0.5
#' 
#' @param max_iter Maximum number of iterations to perform.
#' 
#' @return An R6 object with pip, credible sets, and other features of the fine-mapping result. 
#' \item{pip}{Vector of posterior inclusion probabilities.}
#' 
#' \item{pip_config}{Matrix of posterior inclusion probabilities with each column represents
#' the PIP of corresponding causal configuration. Two ancestries, for eg, the columns represent the 
#' PIP of being first ancestry-specific, second ancestry-specific, and shared.}
#' 
#' \item{alpha}{A length L list of length p by number of causal configuration matrix of posterior inclusion probabilites, with each 
#' column represents the PIP of a scenario.}
#' 
#' \item{KL}{Vector of single-effect KL divergences.}
#' 
#' \item{sigma2}{Residual variance.}
#' 
#' \item{V}{A length L list of N ancestry by N ancestry prior variance estimate.}
#' 
#' \item{elbo}{Vector storing the the evidence lower bound, or
#'   \dQuote{ELBO}, achieved at each iteration of the model fitting
#'   algorithm, which attempts to maximize the ELBO.}
#' 
#' \item{cs}{Estimated credible sets.}
#' @examples
#' library(MESuSiE)
#' data(summary_stat_sd_list)
#' data(R_mat_list)
#' fit = meSuSie_core(R_mat_list,summary_stat_sd_list,L = 10)
#' @import R6 nloptr Rcpp RcppArmadillo Matrix progress
#' @export
meSuSie_core<-function(R_mat_list,summary_stat_list,L,residual_variance=NULL,prior_weights=NULL,ancestry_weight=NULL,optim_method ="optim",estimate_residual_variance =F,max_iter =100,cor_method ="min.abs.corr",cor_threshold=0.5){
  
  cat("*************************************************************\n
  Multiple Ancestry Sum of Single Effect Model (MESuSiE)          \n
   Visit http://www.xzlab.org/software.html For Update            \n
            (C) 2022 Boran Gao, Xiang Zhou                        \n
              GNU General Public License                          \n
*************************************************************") 
  
  time_start<-Sys.time()
  cat("\n# Start data processing for sufficient statistics \n")
  meSuSieData_obj<-meSuSieData$new(R_mat_list,summary_stat_list)
  
  n_snp = nrow(summary_stat_list[[1]])
  n_ancestry = length(summary_stat_list)
  
  if(is.null(prior_weights)){
	prior_weights = rep(1/n_snp,n_snp)
   }
   if(is.null(ancestry_weight)){
    base_fac_ratio = 3
	base_fac = 1/Reduce("+",lapply(seq(1,n_ancestry),function(x)choose(n_ancestry,x)*base_fac_ratio^(n_ancestry-x)))
	ancestry_weight = unlist(lapply(seq(1,n_ancestry),function(x)rep(base_fac*base_fac_ratio^(n_ancestry-x),choose(n_ancestry,x))))
   }
   used_weights = kronecker(prior_weights, t(ancestry_weight), FUN = "*")

  if(is.null(residual_variance)){
    residual_variance = rep(1,n_ancestry)
  }
  
  cat("# Create MESuSiE object \n")
  meSuSieObject_obj<-meSuSieObject$new(n_snp,L,n_ancestry,residual_variance,used_weights,optim_method,estimate_residual_variance,max_iter,names(summary_stat_list))
  cat("# Start data analysis \n")
  
  n_iter = 0
  for (iter in 1:max_iter) {
    
    comp_residual<-meSuSieObject_obj$compute_residual(meSuSieData_obj,meSuSieObject_obj)
    
    
    for (l_index in seq(1,L,1)) {
      
      comp_residual = comp_residual + meSuSieObject_obj$Xr[[l_index]]
      
      SER_res<-single_effect_regression(comp_residual,meSuSieData_obj$XtX.diag, meSuSieObject_obj,l_index) 
      
      meSuSieObject_obj$par_update(SER_res,l_index)
      	  
      meSuSieObject_obj$compute_KL(SER_res,meSuSieObject_obj$compute_SER_posterior_loglik(meSuSieData_obj,comp_residual,SER_res$b1b2),l_index)
      
      meSuSieObject_obj$compute_Xr(meSuSieData_obj,SER_res$b1b2$EB1,l_index)
      
      comp_residual = comp_residual - meSuSieObject_obj$Xr[[l_index]]
         
    }

    updated_sigma2 = meSuSieObject_obj$update_residual_variance(meSuSieData_obj,iter)
    
    if((meSuSieObject_obj$ELBO[iter+1] - meSuSieObject_obj$ELBO[iter])<0.001){
      break
    }
    if(meSuSieObject_obj$estimate_residual_variance==TRUE){
      meSuSieObject_obj$sigma2 =  updated_sigma2
    }
    n_iter = n_iter + 1
    #check convergence and update sigma2
  }

  cat("\n# Data analysis is done, and now generates result \n\n")
  ###Use function in Utility to output result
  
  meSuSieObject_obj$get_result(meSuSie_get_cs(meSuSieObject_obj,R_mat_list,cor_method=cor_method,cor_threshold=cor_threshold),meSusie_get_pip_either(meSuSieObject_obj),meSusie_get_pip_config(meSuSieObject_obj))
  meSuSieObject_obj$mesusie_summary(meSuSieData_obj)

  time_end<-Sys.time()
  cat(c("\n# Total time used for the analysis:",paste0(round(as.numeric(difftime(time_end,time_start,units=c("mins"))),2)," mins\n")))
  return(meSuSieObject_obj)
}


#' MESuSiE Object
#'
#' This object is designed to prepare the MESuSiE analysis.
#'
#' @name meSuSieObject
#' @description An R6 object for handling MESuSiE analysis.
#' @field name_config Character vector for names of the possible configurations.
#' @field column_config List of possible configurations for columns.
#' @field alpha List of alpha matrices for each component.
#' @field mu1 List of mu1 matrices for each component.
#' @field mu2 List of mu2 matrices for each component.
#' @field EB1 List of EB1 matrices for each component.
#' @field EB2 List of EB2 matrices for each component.
#' @field Xr List of Xr matrices for each component.
#' @field KL Numeric vector of KL values for each component.
#' @field lbf Numeric vector of lbf values for each component.
#' @field lbf_variable List of lbf values for each variable.
#' @field ELBO Numeric vector of ELBO values for each iteration.
#' @field sigma2 Numeric vector of residual variances.
#' @field V List of V matrices for each component.
#' @field pi Numeric vector of prior weights.
#' @field estimate_prior_method Character string specifying the method to estimate prior.
#' @field estimate_residual_variance Logical, if TRUE, the residual variance is estimated.
#' @field L Integer for number of causal components.
#' @field nSNP Integer for number of SNPs.
#' @field nancestry Integer for number of ancestries.
#' @field cs List of credible sets.
#' @field pip Numeric vector of PIP values.
#' @field pip_config Matrix of PIP values for each configuration.
#'
#' @description Initialize the MESuSiE object.
#' @param p Integer, number of SNPs.
#' @param L Integer, number of causal components.
#' @param N_ancestry Integer, number of ancestries.
#' @param residual_variance Numeric vector of initial residual variances.
#' @param prior_weights Numeric vector of prior weights.
#' @param estimate_prior_method Character string, method to estimate prior.
#' @param estimate_residual_variance Logical, if TRUE, estimate the residual variance.
#' @param max_iter Integer, maximum number of iterations.
#' @param name_vec Character vector, names for configurations.
#' @return A new `meSuSieObject` object.
#'
#' @description Compute the Xb for given b.
#' @param meSuSie_Data List, MESuSiE data object.
#' @param b Matrix, current estimate of b.
#' @return A list of computed Xb matrices.
#'
#' @description Compute the residual.
#' @param meSuSie_Data List, MESuSiE data object.
#' @param meSuSie_Obj List, current MESuSiE object.
#' @return Matrix of computed residuals.
#'
#' @description Update parameters from the SER results.
#' @param SER_res List, results from the SER analysis.
#' @param l_index Integer, index for the current component.
#' @return Updated parameters in the `meSuSieObject` object.
#' @import R6
#' @export

meSuSieObject <- R6::R6Class("meSuSieObject",public = list(
  initialize = function(p,L,N_ancestry, residual_variance,prior_weights,estimate_prior_method,estimate_residual_variance,max_iter,name_vec){
    
	self$name_config = Reduce(append , lapply(seq(1,length(name_vec)),function(x){
		poss_config = combn(name_vec,x)
		apply(poss_config,2,function(x)paste(x, collapse = '_'))
	}))

	self$column_config = Reduce(append,lapply(seq(1,N_ancestry),function(x){
		poss_config = combn(1:N_ancestry,x)
		lapply(1:ncol(poss_config), function(y)return(matrix(poss_config[,y])))
	}))
	
	self$alpha =rep(list(matrix(0,nrow = p,ncol = N_ancestry)),L)
    self$mu1 = rep(list(lapply(self$column_config,function(x){matrix(0,nrow = p,ncol = length(x))})),L)
    self$mu2 = rep(list(lapply(self$column_config,function(x){matrix(0,nrow = p,ncol = length(x))})),L)
	self$EB1 =rep(list(matrix(0,nrow = p,ncol = N_ancestry)),L)
	self$EB2 =rep(list(matrix(0,nrow = p,ncol = N_ancestry)),L)
	
    self$Xr     = rep(list(matrix(0,nrow = p,ncol=N_ancestry)),L)
    self$KL     = rep(as.numeric(NA),L)
    self$lbf    = rep(as.numeric(NA),L)
    self$lbf_variable = vector("list",L)
    self$ELBO = rep(NA,max_iter)
    self$ELBO[1] = -Inf
    self$sigma2 = residual_variance
    
    
    self$V      = rep(list(matrix(0,ncol = N_ancestry,nrow=N_ancestry)),L)
    self$pi     = prior_weights
    
    self$estimate_prior_method = estimate_prior_method
    self$estimate_residual_variance = estimate_residual_variance
    
    self$L = L
    self$nSNP = p
    self$nancestry = N_ancestry
    
    self$cs = list()
    self$pip = rep(as.numeric(NA),p)
	self$pip_config = matrix(as.numeric(NA),p,2^N_ancestry-1)
  },    
  
  compute_Xb = function (meSuSie_Data,b){
    lapply(1:length(meSuSie_Data$XtX_list),function(x){
      meSuSie_Data$XtX_list[[x]]%*%b[x,]
    })},
  
  compute_residual = function (meSuSie_Data,meSuSie_Obj) {
    residual<-Reduce(cbind,meSuSie_Data$Xty_list)-Reduce("+",meSuSie_Obj$Xr)
    return(residual)
  },
    compute_Xr = function(meSuSie_Data,b1b2,l_index){
    self$Xr[[l_index]] = Reduce(cbind,lapply(1:self$nancestry,function(ancestry_index)meSuSie_Data$XtX_list[[ancestry_index]]%*%b1b2[,ancestry_index]))
  },

  par_update = function(SER_res,l_index){
    self$alpha[[l_index]]<-SER_res$alpha
    self$mu1[[l_index]]<-SER_res$mu1_multi
    self$mu2[[l_index]]<-SER_res$mu2_multi
    self$lbf_variable[[l_index]]<-SER_res$lbf_multi
    self$V[[l_index]]<-SER_res$V
	self$EB1[[l_index]] =SER_res$b1b2$EB1
	self$EB2[[l_index]] =SER_res$b1b2$EB2
	
  },

   compute_SER_posterior_loglik = function(meSuSie_Data,comp_residual,b1b2){
    return(Reduce("+",lapply(1:self$nancestry,function(x){
      sum(-0.5/self$sigma2[[x]]*(-2*comp_residual[,x]*b1b2$EB1[,x]+meSuSie_Data$XtX.diag[[x]]*b1b2$EB2[,x]))
    })))
    
  },
  
  compute_KL = function(SER_res,value,l_index){
    
    self$KL[l_index] = -SER_res$loglik+value
    
  },
  
  update_residual_variance = function(meSuSie_Data,niter){
    

    B_1_ancestry = lapply(1:self$nancestry,function(y)Reduce(cbind,lapply(self$EB1,function(x)x[,y])))
    
    
    BXXB = lapply(1:self$nancestry,function(x){
      sum((t(B_1_ancestry[[x]])%*% meSuSie_Data$XtX_list[[x]])*t(B_1_ancestry[[x]]))
    })  ####{E(BX)}^2
    
    betabar = Reduce("+",self$EB1)
    
    BbarXXBbar = lapply(1:self$nancestry,function(x){
      sum(betabar[,x]*(meSuSie_Data$XtX_list[[x]]%*%betabar[,x]))
    }) ###{E(BbarX)}^2
    
    BbarXty = lapply(1:self$nancestry,function(x){
      2*sum(betabar[,x]*meSuSie_Data$Xty_list[[x]])
    }) ###{E(BbarX)}^2
    

    
    B_2_ancestry = lapply(1:self$nancestry,function(y)Reduce(cbind,lapply(self$EB2,function(x)x[,y])))
    
    XXB2 = lapply(1:self$nancestry,function(x){
      sum(meSuSie_Data$XtX.diag[[x]]*B_2_ancestry[[x]])
    })
    
    updated_sigma = lapply(1:self$nancestry,function(x){
      
      ((meSuSie_Data$yty_list[[x]]-BbarXty[[x]]+BbarXXBbar[[x]])+(XXB2[[x]]-BXXB[[x]]))/meSuSie_Data$N_list[[x]]
    })
    
    
    self$ELBO[niter+1] = Reduce(sum,lapply(1:self$nancestry,function(x){
      -meSuSie_Data$N_list[[x]]/2*log(2*pi*self$sigma2[x])-(1/(2*self$sigma2[x]))*((meSuSie_Data$yty_list[[x]]-BbarXty[[x]]+BbarXXBbar[[x]])+(XXB2[[x]]-BXXB[[x]]))
    }))-sum(self$KL)
    
    return( unlist(updated_sigma))
  },
  
  get_result = function(cs, pip, pip_config){
    self$cs = cs
    self$pip = pip
	self$pip_config = pip_config
	colnames(self$pip_config) = self$name_config
  },
  mesusie_summary = function(meSuSie_Data){
    
    cat(c(paste0("Potential causal SNPs with PIP > 0.5: "),meSuSie_Data$Summary_Stat[[1]]$SNP[which(self$pip>0.5)],"\n\n"))
    cat("Credible sets for effects: \n")
    print(self$cs)
    cat("\n Use MESuSiE_Plot() for visualization")

  }
  
),lock_objects = F
)


##############################################################
#  Assume var(y) = 1 and causal snps contribute negligible variance
#  therefore se(\beta) = var(y)*diag(xtx) => diag(xtx) = var(y)/se(beta)
#                          or
#  R^2 = z^2/(z^2+N-2) => sigma^2 = var(y)*(N-1)/(z^2+N-2)=>diag(xtx)=sigma^2/se(beta)^2
#
#  xtx = sqrt(diag(xtx))%*%R%*%sqrt(diag(xtx))
#  xty = diag(xtx)\beta
#
#############################################################
#' meSuSieData R6 Class
#' 
#' This R6 class is designed to handle and process data relevant to the meSuSie approach.
#' The class includes methods for initialization and processing XtX and XtY data structures.
#' 
#' @field R A list of data.
#' @field Summary_Stat Summary statistics data.
#' @field var_y Variance of y. Default is 1.
#' @field Name_list List of names extracted from `R`.
#' @field N_ancestry Number of ancestry elements in `X`.
#' @field XtX.diag Diagonal elements of XtX.
#' @field XtX_list Processed XtX list.
#' @field Xty_list Processed XtY list.
#' @field N_list List of median values of N from `Summary_Stat`.
#' @field yty_list List of computed yty values.
#' 
#' @description The meSuSieData class is used for handling and processing data in the context
#' of the meSuSie method.
#' @import R6 
#' @export
meSuSieData <- R6::R6Class("meSuSieData",public = list(
  initialize = function(X,Y,var_y = 1){
    self$R<- X
    self$Summary_Stat <-Y
    self$var_y = var_y
    
    
    self$Name_list<-as.list(names(self$R))
    names(self$Name_list)<-names(self$R)
    self$N_ancestry<-length(X)
    
    self$XtX.diag<-self$XtX_diag(self$Summary_Stat,self$Name_list)
    self$XtX_list<-self$XtX_pro(self$R,self$XtX.diag,self$Name_list)
    self$Xty_list<-self$Xty_pro(self$Summary_Stat,self$XtX.diag,self$Name_list) ##diag(xtx)^*betahat
    
    self$N_list<-lapply(self$Summary_Stat,function(x)median(x$N))
    self$yty_list<-lapply(self$N_list,function(x)return(self$var_y*(x-1)))
    
    return(self)},
  ##First compute XtX.diag
  XtX_diag = function(Summary_Stat,Name_list){
    return( lapply(Name_list,function(x){
      R2 = (Summary_Stat[[x]]$Z^2)/(Summary_Stat[[x]]$Z^2+Summary_Stat[[x]]$N-2)
      sigma2 = self$var_y*(1-R2)*(Summary_Stat[[x]]$N-1)/(Summary_Stat[[x]]$N-2)
      return(sigma2/(Summary_Stat[[x]]$Se)^2)
    }))
  },
  
  
  ###Process XtX
  XtX_pro = function(R,XtX.diag,Name_list){
    return(lapply(Name_list,function(x){
      return(diag(sqrt(XtX.diag[[x]]))%*%R[[x]]%*%diag(sqrt(XtX.diag[[x]])))
    }))
  },
  
  ###Process XtY
  Xty_pro = function(Summary_Stat,XtX.diag,Name_list){
    return(lapply(Name_list,function(x){
      XtX.diag[[x]]*Summary_Stat[[x]]$Beta
    }))
    
  }
),
lock_objects = F)

single_effect_regression<-function(XtR,XtX.diag, meSuSieObject_obj,l_index){
  column_config = meSuSieObject_obj$column_config
  N_ancestry = meSuSieObject_obj$nancestry
  Xty_standardized =Reduce(cbind, lapply(1:N_ancestry,function(x){
    XtR[,x]/meSuSieObject_obj$sigma2[x]
  }))
  
  shat2 = Reduce(cbind, lapply(1:N_ancestry,function(x){
    meSuSieObject_obj$sigma2[x]/XtX.diag[[x]]
  }))
  
  betahat = shat2 * Xty_standardized
  
  if(meSuSieObject_obj$estimate_prior_method =="optim"){    
	opt_par<-pre_optim(N_ancestry,-30,10)		
	if(N_ancestry==2){		 
		update_V<-optim(opt_par$inital_par,fn = loglik_cpp,gr=NULL,betahat=betahat,shat2=shat2,prior_weight=meSuSieObject_obj$pi,nancestry =opt_par$nancestry,diag_index = opt_par$diag_index,config_list =column_config,method = "L-BFGS-B",lower=opt_par$lower_bound,upper=opt_par$upper_bound)
		V_mat = vec_to_cov(update_V$par,opt_par$diag_index, opt_par$nancestry)	 
	}else{
		intermediate_V<-nloptr(opt_par$inital_par, eval_f=loglik_cpp,eval_grad_f = NULL,lb = opt_par$lower_bound,ub = opt_par$upper_bound,betahat=betahat,shat2=shat2,prior_weight=meSuSieObject_obj$pi,nancestry =opt_par$nancestry,diag_index = opt_par$diag_index,config_list =column_config,opts=list("algorithm"= "NLOPT_GN_DIRECT_L","xtol_rel"=1.0e-10))
		update_V<-nloptr(intermediate_V$solution, eval_f=loglik_cpp,eval_grad_f = NULL,lb = opt_par$lower_bound,ub = opt_par$upper_bound,betahat=betahat,shat2=shat2,prior_weight=meSuSieObject_obj$pi,nancestry =opt_par$nancestry,diag_index = opt_par$diag_index,config_list=column_config,opts=list("algorithm"= "NLOPT_LN_BOBYQA","xtol_rel"=1.0e-10))
		V_mat = vec_to_cov(update_V$solution,opt_par$diag_index, opt_par$nancestry) 
	}
}

  reg_out<-lapply(column_config,function(x){
	if(length(x)==1){
		uni_reg(betahat[,x],shat2[,x],V_mat[x,x])
	}else if(length(x)>1){
		mvlmm_reg(betahat[,x],shat2[,x],V_mat[x,x]) 
		}
	})

	lbf<-Reduce(cbind,lapply(reg_out,function(x)x$lbf))
	lbf[is.na(lbf)]<-0
	softmax_out<-compute_softmax(lbf,meSuSieObject_obj$pi)
	mu1_multi = lapply(reg_out,function(x)x$post_mean)
	mu2_multi = lapply(reg_out,function(x)x$post_mean2)
	b1b2 = compute_b1b2(softmax_out$alpha_wmulti,mu1_multi,mu2_multi,column_config,ncol(betahat),nrow(betahat))
	
  return(list(alpha = softmax_out$alpha_wmulti,mu1_multi = mu1_multi,mu2_multi = mu2_multi,lbf_multi = lbf,V = V_mat,loglik =softmax_out$loglik,b1b2 = b1b2))
}



