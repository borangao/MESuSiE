
######################################
#
#
#       meSuSie utility functions
#
#
#####################################

#######################################################
#
#
#      Set of functions for single effect model
#
#
#######################################################


#####################################################################################################
# pre_optim is a function to prepare the input for parameter space of the prior covariance matrix
#
# Function takes input of number of ancestry, lower and upper bound of the prior variance
# 
#####################################################################################################
pre_optim<-function(nancestry,min_var,max_var){
  total_ele_num<-(nancestry+1)*nancestry/2
  diag_ele_index<-do.call(rbind,lapply(seq(1,nancestry),function(x)x*(x+1)/2))
  inital_par<-rep(0,total_ele_num)
  lower_bound<-rep(-0.9999,total_ele_num)
  upper_bound<-rep(0.9999,total_ele_num)
  lower_bound[diag_ele_index]<-min_var
  upper_bound[diag_ele_index]<-max_var
  return(list(inital_par=inital_par,lower_bound=lower_bound,upper_bound=upper_bound,diag_index = diag_ele_index,nancestry=nancestry))
}
#####################################################################################################
# vec_to_cov is a function to transform the estimated paramter to covariance matrix
#
# Function takes input of estimated parameter vector, the index of variance and number of ancestry
# 
#####################################################################################################
vec_to_cov<-function(v_vec,diag_index,nancestry){
  cor_mat = matrix(NA,ncol=nancestry,nrow = nancestry)
  cor_mat[upper.tri(cor_mat,diag=TRUE)]<-v_vec
  cor_mat[lower.tri(cor_mat)] <- t(cor_mat)[lower.tri(cor_mat)]
  diag(cor_mat)<-1
  se_mat<-diag(sqrt(exp(v_vec[diag_index])))
  V_mat<-se_mat%*%cor_mat%*%se_mat
 # print(V_mat)
  return(V_mat)
  
}

#####################################################################################################
# multivariate_regression is a function to compute the lbf, posterior mean and posterior second moment
#
# Function takes input of Xty, shat and estimated covariance matrix
# 
#####################################################################################################
uni_reg<-function(betahat,shat2,sigma2){
		lbf = dnorm(betahat,0,sqrt(sigma2 + shat2),log = TRUE) - dnorm(betahat,0,sqrt(shat2),log = TRUE)
		post_var = (1/sigma2 + 1/shat2)^(-1)
		post_mean = post_var * betahat/shat2
		post_mean2 = post_var + post_mean^2
		return(list("lbf" = lbf , post_mean=post_mean,post_mean2=post_mean2))
}

compute_b1b2<-function(alpha,mu1,mu2,column_config,nancestry,nsnp){
	
		EB1 = matrix(0,ncol=nancestry,nrow=nsnp)
		EB2 = matrix(0,ncol=nancestry,nrow=nsnp)
		
		for(x in 1:length(column_config)){
		cor_cols = column_config[[x]]
		EB1[,cor_cols]	= EB1[,cor_cols]+alpha[,x]*mu1[[x]]
		EB2[,cor_cols] = EB2[,cor_cols]+alpha[,x]*mu2[[x]]
		}


	return(list("EB1"=EB1,"EB2"=EB2))

}
multivariate_regression<-function(Xty,shat2,V_mat){
  
  out<-lapply(seq(1,nrow(shat2)),function(x){
    SIGMA<-diag(1/shat2[x,])
    V_mat_SIGMA<-V_mat%*%SIGMA
    post_var<-diag(shat2[x,])-solve(SIGMA+SIGMA%*%V_mat_SIGMA)
    Xty_var <- Xty[x,]%*%post_var	
    post_mean_wmulti<-Xty_var ###Connection of two distribution
    post_mean2_wmulti<-(Xty_var)^2+ diag(post_var)
    
    lbf_part_1<-0.5*sum(Xty_var*Xty[x,])
    lbf_part_2<-0.5*log(det(V_mat_SIGMA+diag(rep(1,ncol(shat2)))))
    return(list(lbf = lbf_part_1 - lbf_part_2,post_mean_wmulti=post_mean_wmulti,post_mean2_wmulti=post_mean2_wmulti))
  })
  return(out)
  
}

#####################################################################################################
# compute_softmax is a function to compute the posterior inclusion probability for l'th effect by lbf
#
# Function takes input of lbf and prior_weights
# 
#####################################################################################################

compute_softmax<-function(lbf,prior_weights){
  maxlbf = max(lbf)
  w_multi = exp(lbf - maxlbf)
  w_weighted_multi = w_multi * prior_weights
  weighted_sum_w = sum(w_weighted_multi)
  alpha_wmulti = w_weighted_multi / weighted_sum_w
  
  return(list(alpha_wmulti = alpha_wmulti, loglik = maxlbf+log(weighted_sum_w)))
}


#######################################################
#
#
#      Set of functions for processing result
#
#
#######################################################

n_in_CS_x = function (x, coverage = 0.9)
  sum(cumsum(sort(x,decreasing = TRUE)) < coverage) + 1
in_CS = function (res, coverage = 0.9) {
  if (inherits(res,"susie"))
    res = res$alpha
  return(t(apply(res,1,function(x) in_CS_x(x,coverage))))
}
in_CS_x = function (x, coverage = 0.9) {
  n = n_in_CS_x(x,coverage)
  o = order(x,decreasing = TRUE)
  result = rep(0,length(x))
  result[o[1:n]] = 1
  return(result)
}
get_purity = function(pos, Xcorr) {
  if (length(pos) == 1)
    c(1,1,1)
  else{
    value = Reduce(cbind,lapply(Xcorr,function(x)c(abs(x[pos,pos]))))
    value = do.call(pmax,data.frame(value))
    return(c(min(value,na.rm = TRUE),
             mean(value,na.rm = TRUE),
             median(value,na.rm = TRUE)))
  }
}
meSuSie_get_cs<-function(res,Xcorr,coverage = 0.95,prior_tol = 1e-9,cor_method =cor_method,cor_threshold=cor_threshold){
  
  include_idx = unlist(lapply(res$V,function(x)max(diag(x))>prior_tol))
  alpha_either = t(Reduce(cbind,lapply(res$alpha,function(x)apply(x,1,sum))))
  status = in_CS(alpha_either,coverage = 0.95)
  cs = lapply(1:nrow(status),function(i) which(status[i,]!=0))
  claimed_coverage = sapply(1:length(cs),function (i) sum(alpha_either[i,][cs[[i]]]))
  include_idx = include_idx * (lapply(cs,length) > 0)
  include_idx = include_idx * (!(duplicated(cs)))
  
  include_idx = as.logical(include_idx)
  if (sum(include_idx) == 0)
    return(list(cs = NULL,coverage = NULL,requested_coverage = coverage))
  
  cs = cs[include_idx]
  claimed_coverage = claimed_coverage[include_idx]
  
  purity = data.frame(do.call(rbind,lapply(1:length(cs),function (i){
    get_purity(cs[[i]],Xcorr)
  })))
  
  colnames(purity) = c("min.abs.corr","mean.abs.corr","median.abs.corr")
  threshold = cor_threshold
  is_pure = which(purity[,colnames(purity)%in%cor_method] >= threshold)
  
 # is_pure = which(purity[,1] >= threshold)
  
  if (length(is_pure) > 0) {
    cs = cs[is_pure]
    purity = purity[is_pure,]
    row_names = paste0("L",which(include_idx)[is_pure])
    names(cs) = row_names
    rownames(purity) = row_names
    
    # Re-order CS list and purity rows based on purity.
    ordering = order(purity[,1],decreasing = TRUE)
    return(list(cs       = cs[ordering],
                purity   = purity[ordering,],
                cs_index = which(include_idx)[is_pure[ordering]],
                coverage = claimed_coverage[ordering],
                requested_coverage=coverage))
  } else
    return(list(cs = NULL,coverage = NULL, requested_coverage = coverage))
}

meSusie_get_pip_either<-function(res){
	return(do.call(pmax,data.frame(Reduce(cbind,lapply(res$alpha,function(x)apply(x,1,sum))))))
}
meSusie_get_pip_config = function(res) {
	return(Reduce(cbind,lapply(1:ncol(res$alpha[[1]]),function(x){
		do.call(pmax,lapply(res$alpha,function(y)y[,x]))
	})))
}
meSusie_get_pip = function (res, prior_tol = 1e-9) {
  include_idx = unlist(lapply(res$V,function(x)max(diag(x))>prior_tol))
  
  alpha_either = t(Reduce(cbind,lapply(res$alpha,function(x)apply(x,1,sum))))
  
  
  if (length(include_idx) > 0)
    alpha_include = alpha_either[include_idx,,drop = FALSE]
  else
    alpha_include = matrix(0,1,ncol(res$alpha))
  return(as.vector(1 - apply(1 - alpha_include,2,prod)))
}
######################################
#
#
#       meSuSie Plot Part
#
#
#####################################
require(ggplot2)
require(cowplot)
meSusie_plot_pip<-function(res,R_mat,summary_data,annotation_data=NULL){
  
  
  
  names(summary_data)
#  summary_data<-summary_stat_list
 # summary_data<-lapply(summary_data,function(x){
 #   x$POS=seq(1,500,1)
 #   return(x)})
  ##Create Mahattan Plot for GWAS Data
  
  for(x in 1:res$nancestry){
    data_plot = summary_data[[x]]
    data_plot$Pvalue = 2*pnorm(-abs(data_plot$Z))
    Lead_SNP = data_plot$SNP[which.max(abs(data_plot$Z))]
    Lead_SNP_r2 = R_mat[[x]][,which(colnames(R_mat[[x]])==Lead_SNP)]
    p_manhattan = ggplot() + geom_point(data =data_plot[-(data_plot$SNP==Lead_SNP),], aes(x = POS, y = -log10(Pvalue), color = Lead_SNP_r2[-(data_plot$SNP==Lead_SNP)]),size = 1.5)+
      geom_point(data =data_plot[data_plot$SNP==Lead_SNP,], aes(x = POS, y = -log10(Pvalue)),shape =3,size = 3,color="red")+
      geom_text(data =data_plot[data_plot$SNP==Lead_SNP,], mapping=aes(x=POS, y=-log10(Pvalue), label=SNP),vjust=1.2, size=4,show.legend = FALSE)+
      scale_color_stepsn(colors = c("navy", "lightskyblue", "green", "orange", "red"),breaks = seq(0.2, 0.8, by = 0.2),limits = c(0, 1),show.limits = TRUE,na.value = 'grey50')+
      geom_hline(yintercept = -log10(5e-8),linetype = "dashed",color = "grey50",size = 0.5)+
      geom_vline(xintercept = data_plot$POS[data_plot$SNP==Lead_SNP],linetype = "dashed",color = "grey50",size = 0.5)+
      xlab(paste0(names(summary_data)[x]," GWAS"))+ylab("-log(Pvalue)")+
      theme_bw()+ theme(legend.position="none") #guides(fill=guide_legend(title="r2"))+labs(colour= "r2")
    assign(paste0("GWAS_ancestry_",x),p_manhattan)
  }
  
  SNP_pip_05 = which(res$pip>0.5)
  all_SNP_cs = unlist(res$cs$cs)
  following_snp = setdiff(all_SNP_cs,SNP_pip_05)
  
  SNP_pip_05_label = rep(0,length(res$pip))
  SNP_pip_05_label[SNP_pip_05]=1
  SNP_following_snp_label = rep(0,length(res$pip))
  SNP_following_snp_label[following_snp]=1
  
  data_plot = summary_data[[1]]
  data_plot$Pvalue = res$pip
  data_plot$LEAD_SNP_LABEL = SNP_pip_05_label
  data_plot$following_snp_label =  SNP_following_snp_label
  
  p_pip = ggplot() + geom_point(data =data_plot[(data_plot$LEAD_SNP_LABEL!=1)&(data_plot$following_snp_label!=1),], aes(x = POS, y = Pvalue),size = 1.5,color="gray")+
    geom_point(data =data_plot[data_plot$LEAD_SNP_LABEL==1,], aes(x = POS, y = Pvalue),shape =3,size = 3,color="red")+
    geom_text(data =data_plot[data_plot$LEAD_SNP_LABEL==1,], mapping=aes(x=POS, y=Pvalue, label=SNP),vjust=1.2, size=4,show.legend = FALSE)+
    geom_point(data =data_plot[data_plot$following_snp_label==1,], aes(x = POS, y = Pvalue),size = 1,color="purple")+
    geom_hline(yintercept = 0.5,linetype = "dashed",color = "grey50",size = 0.5)+
    xlab(paste0("meSuSie Result"))+ylab("PIP")+
    theme_bw()
  
  plotlist<-list()
  for(i in 1:res$nancestry){
    
    plotlist[[i]]=get(paste0("GWAS_ancestry_",i))
    
  }
  plotlist[[res$nancestry+1]]<-p_pip
  cowplot::plot_grid(plotlist=plotlist)
  
}


