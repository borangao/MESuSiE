
######################################
#
#
#       meSuSie utility functions
#
#
#####################################

#######################################################################################################################################
#
#
#      Set of functions for single effect model
#
#
#######################################################################################################################################


##########################################################################################
# Function: pre_optim
#
# Purpose: 
# The function prepares the input for the parameter space of the prior covariance matrix.
#
# Parameters:
# - nancestry: Number of ancestries.
# - min_var: Lower bound of the prior variance.
# - max_var: Upper bound of the prior variance.
#
# Returns:
# A list containing the initial parameters, lower bounds, upper bounds, diagonal element 
# indices, and the number of ancestries.
##########################################################################################

pre_optim <- function(nancestry, min_var, max_var) {
  
  # Compute the total number of elements in the covariance matrix
  total_ele_num <- (nancestry + 1) * nancestry / 2
  
  # Compute the indices of diagonal elements in the covariance matrix
  diag_ele_index <- do.call(rbind, lapply(seq(1, nancestry), function(x) x * (x + 1) / 2))
  
  # Set initial parameters to zeros
  inital_par <- rep(0, total_ele_num)
  
  # Set default bounds for off-diagonal elements
  lower_bound <- rep(-0.9999, total_ele_num)
  upper_bound <- rep(0.9999, total_ele_num)
  
  # Override bounds for diagonal elements using provided min_var and max_var
  lower_bound[diag_ele_index] <- min_var
  upper_bound[diag_ele_index] <- max_var
  
  # Return a list containing the constructed parameters and bounds
  return(list(
    inital_par = inital_par,
    lower_bound = lower_bound,
    upper_bound = upper_bound,
    diag_index = diag_ele_index,
    nancestry = nancestry
  ))
}

##########################################################################################
# Function: vec_to_cov
#
# Purpose: 
# Transforms the estimated parameter vector into a covariance matrix.
#
# Parameters:
# - v_vec: The estimated parameter vector.
# - diag_index: Index of the variance in the parameter vector.
# - nancestry: Number of ancestries.
#
# Returns:
# The constructed covariance matrix.
##########################################################################################

vec_to_cov <- function(v_vec, diag_index, nancestry) {
  
  # Initialize an empty correlation matrix filled with NAs
  cor_mat <- matrix(NA, ncol = nancestry, nrow = nancestry)
  
  # Fill the upper triangle (including diagonal) of the correlation matrix using the parameter vector
  cor_mat[upper.tri(cor_mat, diag = TRUE)] <- v_vec
  
  # Mirror the upper triangle values to the lower triangle to make the matrix symmetric
  cor_mat[lower.tri(cor_mat)] <- t(cor_mat)[lower.tri(cor_mat)]
  
  # Set the diagonal of the correlation matrix to 1
  diag(cor_mat) <- 1
  
  # Compute the standard error matrix using the diagonal variance values
  se_mat <- diag(sqrt(exp(v_vec[diag_index])))
  
  # Compute the covariance matrix by combining the correlation and standard error matrices
  V_mat <- se_mat %*% cor_mat %*% se_mat
  
  return(V_mat)
}


##########################################################################################
# Function: uni_reg
#
# Purpose: 
# Computes the log Bayes Factor (lbf) and properties of the posterior distribution for
# given input parameters.
#
# Parameters:
# - betahat: Estimated regression coefficient.
# - shat2: Estimated variance of the coefficient.
# - sigma2: Prespecified variance value.
#
# Returns:
# A list containing the log Bayes Factor (lbf), the mean of the posterior distribution 
# (post_mean), and the second moment of the posterior distribution (post_mean2).
##########################################################################################

uni_reg <- function(betahat, shat2, sigma2) {
  
  # Compute the log Bayes Factor (lbf)
  lbf <- dnorm(betahat, 0, sqrt(sigma2 + shat2), log = TRUE) - 
         dnorm(betahat, 0, sqrt(shat2), log = TRUE)
  
  # Compute the posterior variance
  post_var <- 1 / (1/sigma2 + 1/shat2)
  
  # Compute the mean of the posterior distribution
  post_mean <- post_var * betahat / shat2
  
  # Compute the second moment of the posterior distribution
  post_mean2 <- post_var + post_mean^2
  
  # Return a list containing computed values
  return(list(lbf = lbf, post_mean = post_mean, post_mean2 = post_mean2))
}
##########################################################################################
# Function: compute_b1b2
#
# Purpose: 
# Computes the EB1 and EB2 matrices based on the provided input parameters.
#
# Parameters:
# - alpha: Input alpha matrix.
# - mu1: Input mu1 list.
# - mu2: Input mu2 list.
# - column_config: Configuration of columns.
# - nancestry: Number of ancestries.
# - nsnp: Number of SNPs.
#
# Returns:
# A list containing the computed EB1 and EB2 matrices.
##########################################################################################

compute_b1b2<-function(alpha,mu1,mu2,column_config,nancestry,nsnp){
		# Initialize matrices with zeros
		EB1 = matrix(0,ncol=nancestry,nrow=nsnp)
		EB2 = matrix(0,ncol=nancestry,nrow=nsnp)
		# Loop over column configurations
		for(x in 1:length(column_config)){
		cor_cols = column_config[[x]]
		# Update EB1 and EB2 matrices based on column configurations
		EB1[,cor_cols]	= EB1[,cor_cols]+alpha[,x]*mu1[[x]]
		EB2[,cor_cols] = EB2[,cor_cols]+alpha[,x]*mu2[[x]]
		}

	return(list(EB1 = EB1, EB2 = EB2))

}

##########################################################################################
# Function: compute_softmax
#
# Purpose: 
# Computes the posterior inclusion probability for the l'th effect using log Bayes Factors.
#
# Parameters:
# - lbf: Log Bayes Factors.
# - prior_weights: A vector of prior weights.
#
# Returns:
# A list containing the posterior inclusion probability (alpha_wmulti) and 
# the logarithm of the likelihood (loglik).
##########################################################################################

compute_softmax <- function(lbf, prior_weights) {
  
  # Subtracting max lbf from all lbf values for numerical stability
  maxlbf = max(lbf)
  w_multi = exp(lbf - maxlbf)
  
  # Multiply by prior weights
  w_weighted_multi = w_multi * prior_weights
  
  # Sum over weighted values
  weighted_sum_w = sum(w_weighted_multi)
  
  # Compute the normalized weights
  alpha_wmulti = w_weighted_multi / weighted_sum_w
  
  # Return computed values
  return(list(alpha_wmulti = alpha_wmulti, loglik = maxlbf + log(weighted_sum_w)))
}


#######################################################################################################################################
#
#
#      Set of functions for processing result
#
#
#######################################################################################################################################

##########################################################################################
# Function: n_in_CS_x
#
# Purpose: 
# Calculates the number of SNPs in a credible set required to achieve the desired coverage.
#
# Parameters:
# - x: A numeric vector representing SNP values.
# - coverage: The desired coverage (default is 0.95 or 95%).
#
# Returns:
# The number of SNPs needed in the credible set.
##########################################################################################

n_in_CS_x <- function (x, coverage = 0.95) {
  
  # Cumulatively sum the sorted SNPs in descending order 
  # and find the position just before the cumulative sum exceeds the desired coverage
  num_snps_required = sum(cumsum(sort(x, decreasing = TRUE)) < coverage) + 1
  
  return(num_snps_required)
}

##########################################################################################
# Function: in_CS_x
#
# Purpose: 
# Identifies the SNPs that are included in the credible set to achieve the desired coverage.
#
# Parameters:
# - x: A numeric vector representing SNP values.
# - coverage: The desired coverage (default is 0.95 or 95%).
#
# Returns:
# A binary vector, where 1 indicates the SNP is in the credible set and 0 otherwise.
##########################################################################################

in_CS_x <- function (x, coverage = 0.95) {
  
  # Get the number of SNPs required for the desired coverage
  n = n_in_CS_x(x, coverage)
  
  # Order the SNPs in descending order
  o = order(x, decreasing = TRUE)
  
  # Initialize result vector as zeros
  result = rep(0, length(x))
  
  # Mark the top 'n' SNPs as 1, indicating they are in the credible set
  result[o[1:n]] = 1
  
  return(result)
}

##########################################################################################
# Function: in_CS
#
# Purpose: 
# Applies the in_CS_x function row-wise to a matrix to identify the SNPs that are 
# included in the credible set for each row.
#
# Parameters:
# - res: A matrix where each row represents a set of SNP values.
# - coverage: The desired coverage (default is 0.95 or 95%).
#
# Returns:
# A matrix where each row is a binary vector indicating which SNPs are in the credible set.
##########################################################################################

in_CS <- function (res, coverage = 0.95) {
  
  # Apply the in_CS_x function to each row of the 'res' matrix
  result_matrix = t(apply(res, 1, function(x) in_CS_x(x, coverage)))
  
  return(result_matrix)
}

##########################################################################################
# Function: get_purity
#
# Purpose: 
# Computes the minimum, mean, and median of absolute correlations for specific positions.
#
# Parameters:
# - pos: A vector of positions to consider.
# - Xcorr: A correlation matrix.
#
# Returns:
# A vector containing the minimum, mean, and median of the absolute correlations.
##########################################################################################

get_purity <- function(pos, Xcorr) {
  
  # If there's only one position, return a vector of ones
  if (length(pos) == 1) {
    return(c(1, 1, 1))
  } else {
    # Extract absolute correlations for specified positions
    value_list = lapply(Xcorr, function(x) c(abs(x[pos, pos])))
    
    # Combine the lists column-wise
    value_matrix = Reduce(cbind, value_list)
    
    # Compute the maximum for each row
    value_max = do.call(pmax, data.frame(value_matrix))
    
    # Return the desired statistics
    return(c(min(value_max, na.rm = TRUE),
             mean(value_max, na.rm = TRUE),
             median(value_max, na.rm = TRUE)))
  }
}
#' Extract Credible Sets from MESuSiE Results
#'
#' This function processes results from MESuSiE to identify credible sets,
#' compute their purity and provide other associated statistics.
#'
#' @param res The result object from MESuSiE.
#' @param Xcorr A list of correlation matrix used in MESuSiE.
#' @param coverage A numeric value indicating the desired coverage for the credible set, default is \(0.95\).
#' @param prior_tol A numeric threshold value to determine prior variances, default is \(1e-9\).
#' @param cor_method A string indicating which method to use when computing purity, defalt is "min.abs.corr".
#' @param cor_threshold A numeric threshold value for purity measurement, default is 0.5.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{cs} The credible sets.
#'   \item \code{cs_category} The category of each credible set being shared or ancestry specific.
#'   \item \code{purity} Purity statistics for each credible set.
#'   \item \code{cs_index} The index for the credible sets.
#'   \item \code{coverage} The achieved coverage for the credible sets.
#'   \item \code{requested_coverage} The desired coverage (provided as input).
#' }
#' @export

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
    return(list(cs = NULL,cs_category = NULL,purity = NULL, cs_index = NULL,coverage = NULL, requested_coverage = coverage))
  
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
    claimed_coverage = claimed_coverage[is_pure]
    row_names = paste0("L",which(include_idx)[is_pure])
    names(cs) = row_names
    rownames(purity) = row_names
    
    cs_index = which(include_idx)[is_pure]
    cs_category = unlist(lapply(res$alpha[cs_index],function(x)res$name_config[which.max(colSums(x))]))
    names(cs_category) = row_names
    
    return(list(cs       = cs,
                cs_category = cs_category,
                purity   = purity,
                cs_index = cs_index,
                coverage = claimed_coverage,
                requested_coverage=coverage))
  } else
    return(list(cs = NULL,cs_category = NULL,purity = NULL, cs_index = NULL,coverage = NULL, requested_coverage = coverage))
}

##########################################################################################
# Function: meSusie_get_pip_either
#
# Purpose: 
# Extract the maximum posterior inclusion probabilities (PIP) across all subsets from the given results.
#
# Parameters:
# - res: Results likely containing alpha values representing inclusion probabilities.
#
# Returns:
# A vector containing the maximum PIP for each result.
##########################################################################################

meSusie_get_pip_either <- function(res) {
  
  # Aggregate the alpha values for each subset
  aggregated_alpha = Reduce(cbind, lapply(res$alpha, function(x) apply(x, 1, sum)))
  
  # Return the maximum value for each result
  return(do.call(pmax, data.frame(aggregated_alpha)))
}

##########################################################################################
# Function: meSusie_get_pip_config
#
# Purpose: 
# Extract the maximum posterior inclusion probabilities (PIP) for each configuration from the given results.
#
# Parameters:
# - res: Results containing alpha matrices with inclusion probabilities.
#
# Returns:
# A matrix containing the maximum PIP for each configuration across all results.
##########################################################################################

meSusie_get_pip_config <- function(res) {
  
  # Extract the maximum alpha value for each column across all alpha matrices
  max_alpha_per_config = Reduce(cbind, lapply(1:ncol(res$alpha[[1]]), function(x) {
    do.call(pmax, lapply(res$alpha, function(y) y[,x]))
  }))
  
  return(max_alpha_per_config)
}

#######################################################################################################################################
#
#
#      MESuSiE Plot Part
#
#
#######################################################################################################################################
custom_theme <- function() {
  theme(
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 7),  
    axis.title.x = element_text(size = 7, face="bold"),
    axis.title.y = element_text(size = 7, face="bold"),
    strip.text.x = element_text(size = 7),
    strip.text.y = element_text(size = 7),
    strip.background = element_blank(),
    legend.text = element_text(size=7),
    legend.title = element_text(size=7, face="bold"),
    plot.title = element_text(size=7, hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  )
}
#' Create GWAS and PIP Plots from MESuSiE Results
#'
#' This function generates GWAS (Genome-Wide Association Study) and PIP (Posterior Inclusion Probability)
#' plots for each ancestry based on the results from MESuSiE.
#'
#' @param res The result object from MESuSiE.
#' @param R_mat A list containing r^2 matrices for each ancestry.
#' @param summary_data A list of data frames for each ancestry containing at least SNP, Z (Z-scores), 
#'     and possibly POS (position) columns.
#'
#' @return A combined ggplot object of GWAS plots for each ancestry and a PIP plot.
#'
#' @examples
#' \dontrun{
#' # Assuming `results` is from MESuSiE, `ld_matrices` is the LD matrix list and `data_summary` is the summary data list
#' MESuSiE_Plot_out <- MESuSiE_Plot(results, ld_matrices, data_summary)
#' print(MESuSiE_Plot_out)
#' }
#'
#' @export
#' @importFrom dplyr select pull mutate bind_cols
#' @importFrom tidyr pivot_longer
#' @import ggplot2 
#' @import cowplot 
#' @import ggrepel 
MESuSiE_Plot <- function(res, R_mat, summary_data) {
  
  # Extract names of the datasets in summary_data
  name_vec <- names(summary_data)
  
  ### 1. Position Information Check ###
  
  # Check if Position information is provided; if not, generate a position based on row order
  if (!("POS" %in% colnames(summary_data[[name_vec[1]]]))) {
    summary_data <- lapply(summary_data, function(x) {
      x$POS = seq_len(nrow(x))
      return(x)
    })
  }
  
  # Extract baseline data (SNP and POS columns)
  Baseline_data <- summary_data[[names(summary_data)[1]]] %>% select(SNP, POS)
  
  ### 2. Construct Z_data ###
  
  # Combine Z values across different datasets in summary_data
  Z_data <- Reduce(cbind, lapply(name_vec, function(x) {
    Z <- matrix(summary_data[[x]] %>% pull(Z), ncol = 1)
    colnames(Z) <- paste0("Z_", x)
    return(Z)
  }))
  
  # Extract PIP information
  PIP_data <- data.frame(PIP = res$pip)
  
  # Combine Baseline_data, Z_data, and PIP_data
  all_data <- bind_cols(Baseline_data, Z_data, PIP_data)
  
  ### 3. Identify Lead SNP and Ancestry ###
  
  # Determine which SNP and ancestry has the maximum Z value
  lead_SNP_Ancestry <- all_data %>%
    select(SNP, starts_with("Z_")) %>%
    pivot_longer(cols = -SNP, names_to = "ancestry", values_to = "z_val") %>%
    arrange(desc(abs(z_val))) %>%
    slice(1) %>%
    select(SNP, ancestry)
  
  lead_SNP <- lead_SNP_Ancestry %>% pull(SNP)
  lead_Ancestry <- gsub("Z_", "", lead_SNP_Ancestry %>% pull(ancestry))
  
  ### 4. Calculate LD (r-squared) ###
  
  # Calculate LD r-squared with respect to the lead SNP for all ancestries
  r_data <- Reduce(cbind, lapply(name_vec, function(x) {
    r <- matrix(R_mat[[x]][, lead_SNP], ncol = 1)
    colnames(r) <- paste0("r_", x)
    return(r)
  }))
  
  all_data <- bind_cols(all_data, r_data)
  
  ### 5. GWAS and PIP Plots ###
  
  plotlist <- list()
  
  # Generate GWAS plot for each ancestry
  for (ancestry_name in name_vec) {
    gwas_plot_data <- all_data %>%
      select(SNP, POS, paste0("Z_", ancestry_name), paste0("r_", ancestry_name)) %>%
      mutate(P = -log10(2 * pnorm(-abs(!!sym(paste0("Z_", ancestry_name))))))
    
    locus_zoom <- ggplot(gwas_plot_data, aes(x = SNP, y = P, color = !!sym(paste0("r_", ancestry_name)))) +
      geom_point(size = 1.5) +
      scale_color_stepsn(
        colors = c("navy", "lightskyblue", "green", "orange", "red"),
        breaks = seq(0.2, 0.8, by = 0.2),
        limits = c(0, 1),
        show.limits = TRUE,
        na.value = 'grey50',
        name = expression(R^2)
      ) +
      xlab(ancestry_name) + ylab("-log10 P-value") +
      geom_hline(yintercept = -log10(5e-8), linetype = "dashed") +
      custom_theme()  +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            legend.position = "none")
    
    plotlist[[ancestry_name]] <- locus_zoom
  }
  
  # Generate a separate PIP plot
  PIP_plot_data <- all_data %>% select(SNP, POS, PIP, paste0("r_", lead_Ancestry))
  PIP_Plot <- ggplot(PIP_plot_data, aes(x = SNP, y = PIP, color = !!sym(paste0("r_", lead_Ancestry)))) +
    geom_point(size = 1.5) +
    scale_color_stepsn(
      colors = c("navy", "lightskyblue", "green", "orange", "red"),
      breaks = seq(0.2, 0.8, by = 0.2),
      limits = c(0, 1),
      show.limits = TRUE,
      na.value = 'grey50',
      name = expression(R^2)
    ) +
    xlab("MESuSiE") + ylab("PIP") +
    geom_hline(yintercept = 0.5, linetype = "dashed") +      
	custom_theme()  +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            legend.position = "bottom")
  legend_used <- get_legend(PIP_Plot)
  PIP_Plot <- PIP_Plot + theme(legend.position = "none")
  
  plotlist[["PIP_Plot"]] <- PIP_Plot
  
  # Combine all plots and return
  combined_plot <- cowplot::plot_grid(plotlist = plotlist, ncol = 1)
  combined_plot_out <- plot_grid(combined_plot, legend_used, ncol = 1, rel_heights = c(1, .1))
  
  return(combined_plot_out)
}

