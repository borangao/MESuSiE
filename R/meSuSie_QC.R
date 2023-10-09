############################################################################
#
#		This file defines the QC function for GWAS, Reference SNP alignment, 
#		and LD mismatch function that can be used in the multiple-ancestry
#		fine-mapping
#
############################################################################
##############################################
#	Step 1. QC of GWAS for each ancestry
###############################################

#' GWAS Quality Control Function
#'
#' This function performs several quality control steps on GWAS data:
#' 1. Removes strand ambiguous variants.
#' 2. Removes multi-allelic variants.
#' 3. Excludes the MHC complex region.
#' 4. Filters SNPs based on Minor Allele Frequency (MAF).
#'
#' @param dat A data frame containing the GWAS data. Columns "CHR", "POS", "CHR_POS", 
#' "REF", "ALT", and (optionally) "MAF" are required.
#' @param MAF_threshold The threshold for Minor Allele Frequency filtering. SNPs 
#' with MAF below this threshold or above (1 - threshold) will be removed.
#'
#' @return A filtered data frame after performing the quality control steps.
#' @import dplyr
#' @export

GWAS_QC <- function(dat, MAF_threshold) { 
	# CHR, POS , CHR_POS, REF, ALT, and MAF are required for QC
	# QC Step 1: Remove Strand Ambiguous Variants
	# These are SNPs that can't be differentiated based on the strand and may introduce errors.
	dat <- dat %>%
	filter(!(
	(REF == "G" & ALT == "C") |
	(REF == "C" & ALT == "G") |
	(REF == "T" & ALT == "A") |
	(REF == "A" & ALT == "T")
	))
  
  # QC Step 2: Remove Multi-allelic Variants
  # These are variants with more than two alleles.
	dat <- dat %>% 
	filter(!(nchar(REF) + nchar(ALT) > 2))	
  
  # QC Step 3: Exclude MHC Complex Region
  # This region can be particularly complex and may need special handling.
  dat <- dat %>% 
    filter(!(CHR == 6 & POS > 25e6 & POS < 34e6))
  
  # QC Step 4: Minor Allele Frequency (MAF) Filter
  # Remove SNPs with MAF below the threshold or above (1 - threshold).
  if("MAF"%in%colnames(dat)){
  dat <- dat %>% 
    filter(MAF > MAF_threshold, MAF < 1 - MAF_threshold)
  }
  dat <- dat %>%
  distinct(CHR_POS,.keep_all=TRUE)
  return(dat)
}
##############################################
#	Step 2. Find common set of SNPs, and align 
#	the reference allele, switch the sign of beta
#	or z score based on the allele
###############################################

#' Identify Common SNPs Between Two Datasets
#'
#' This function identifies common SNPs (based on "CHR_POS") between a summary statistics 
#' dataset (`sumstat`) and a reference dataset (`ref`). It then ensures that the 
#' alleles match between the two datasets, either in the original or flipped orientation.
#'
#' @param sumstat A data frame containing the summary statistics data with at least 
#' the columns "CHR_POS", "REF", and "ALT".
#' @param ref A data frame serving as the reference, with at least the columns 
#' "CHR_POS", "REF", and "ALT".
#'
#' @return A vector of "CHR_POS" values corresponding to the SNPs common between 
#' `sumstat` and `ref` with matching allele information.
#' @import dplyr
#' @export

find_common_snps <- function(sumstat, ref) {
  
  # Identifying common SNPs between sumstat and ref datasets
  common_SNP <- intersect(sumstat$CHR_POS, ref$CHR_POS)
  
  # Subsetting data frames based on common SNPs and arranging by order of common_SNP
  sumstat <- sumstat %>% 
    filter(CHR_POS %in% common_SNP) %>% 
    arrange(match(CHR_POS, common_SNP))
  
  ref <- ref %>% 
    filter(CHR_POS %in% common_SNP) %>% 
    arrange(match(CHR_POS, common_SNP))
  
  matched_pos<-which((sumstat$REF==ref$REF&sumstat$ALT==ref$ALT)|(sumstat$REF==ref$ALT&sumstat$ALT==ref$REF))
  
  common_SNP<-sumstat[matched_pos,]%>%pull(CHR_POS)
  
  return(common_SNP)
}
# Match the reference allele with the sumstat used, 
# flip the beta and zscore if the reference allele is flipped.

#' Allele Flip Function for SNP Matching
#'
#' This function modifies the reference dataset (`ref`) based on the summary 
#' statistics dataset (`sumstat`) to ensure that the alleles match. If the alleles 
#' do not match, the function will flip the alleles and corresponding statistics 
#' in the reference dataset.
#'
#' @param sumstat A data frame containing the summary statistics data with at least 
#' the columns "REF".
#' @param ref A data frame serving as the reference with at least the columns 
#' "REF", "ALT", "BETA", and "SE".
#'
#' @return A modified version of `ref` where alleles have been flipped to match 
#' `sumstat` and corresponding statistics (BETA and Z) have been updated accordingly.
#' @import dplyr
#' @export

allele_flip<-function(sumstat,ref){
	ref<-ref %>% mutate(flip_index = ifelse(ref$REF==sumstat$REF,1,-1))%>%
			 mutate(BETA = BETA * flip_index)%>%
			 mutate(Z = BETA/SE, REF = sumstat$REF, ALT = sumstat$ALT)%>%
			 select(-flip_index)
	return(ref)
}
#################################################
#	Step 3. Match up the reference panel with
#	the GWAS reference and alternative
#	alleles
###############################################
# Here we use plink to do the referernce allele flip
# Write out the reference SNP for the genotype data
# plink --bfile ref_geno --ref-allele ref.txt --make-bed --out ref_geno_flip
################################################
#	Step 4. Check the LD mismatch within each 
#	ancestry. We adpoted kring_rss
#	function of SuSiE to detect potential LD
#	mismatch, we slightly modified the function
#	to disclose SNPs that are not marginally 
#	significant but in high LD with GWAS signal
#	and have a high loglik
###############################################

#' LD mismatch detection (Adapted from SuSiE)
#'
#' This function performs LD mismatch detection on summary statistics and LD
#'
#' @param z A numeric vector of summary statistics.
#' @param R A correlation matrix. If the eigen decomposition hasn't been 
#' performed on R, it will be computed within the function.
#' @param r_tol A numeric tolerance level for positive semidefiniteness 
#' of R. Default is \( 1e-08 \).
#' @param s The parameter \( s \) for RSS-based fine-mapping. 
#' If not provided, it will be estimated using `estimate_s_rss` function 
#' with the "null-mle" method.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{plot} A ggplot object visualizing the observed vs. expected 
#'   z scores.
#'   \item \code{conditional_dist} A data frame containing the summary 
#'   statistics, conditional means, variances, standardized differences, 
#'   and log likelihood ratios.
#' }
#' @importFrom mixsqp mixsqp
#' @export

kriging_rss = function (z, R, r_tol = 1e-08, s = estimate_s_rss(z, R, r_tol, 
    method = "null-mle")) 
{
    if (is.null(attr(R, "eigen"))) 
        attr(R, "eigen") = eigen(R, symmetric = TRUE)
    eigenld = attr(R, "eigen")
    if (any(eigenld$values < -r_tol)) 
        warning("The matrix R is not positive semidefinite. Negative ", 
            "eigenvalues are set to zero.")
    eigenld$values[eigenld$values < r_tol] = 0
    if (s > 1) {
        warning("The given s is greater than 1. We replace it with 0.8.")
        s = 0.8
    }
    if (s < 0) 
        stop("The s must be non-negative")
    dinv = 1/((1 - s) * eigenld$values + s)
    dinv[is.infinite(dinv)] = 0
    precision = eigenld$vectors %*% (t(eigenld$vectors) * dinv)
    condmean = rep(0, length(z))
    condvar = rep(0, length(z))
    for (i in 1:length(z)) {
        condmean[i] = -(1/precision[i, i]) * precision[i, -i] %*% 
            z[-i]
        condvar[i] = 1/precision[i, i]
    }
    z_std_diff = (z - condmean)/sqrt(condvar)
    a_min = 0.8
    if (max(z_std_diff^2) < 1) {
        a_max = 2
    }
    else {
        a_max = 2 * sqrt(max(z_std_diff^2))
    }
    npoint = ceiling(log2(a_max/a_min)/log2(1.05))
    a_grid = 1.05^((-npoint):0) * a_max
    sd_mtx = outer(sqrt(condvar), a_grid)
    matrix_llik = dnorm(z - condmean, sd = sd_mtx, log = TRUE)
    lfactors = apply(matrix_llik, 1, max)
    matrix_llik = matrix_llik - lfactors
    w = mixsqp(matrix_llik, log = TRUE, control = list(verbose = FALSE))$x
    logl0mix = as.numeric(log(exp(matrix_llik) %*% w)) + lfactors
    matrix_llik = dnorm(z + condmean, sd = sd_mtx, log = TRUE)
    lfactors = apply(matrix_llik, 1, max)
    matrix_llik = matrix_llik - lfactors
    logl1mix = as.numeric(log(exp(matrix_llik) %*% w)) + lfactors
    logLRmix = logl1mix - logl0mix
    res = data.frame(z = z, condmean = condmean, condvar = condvar, 
        z_std_diff = z_std_diff, logLR = logLRmix)
    p = ggplot(res) + geom_point(aes(y = z, x = condmean)) + 
        labs(y = "Observed z scores", x = "Expected value") + 
        geom_abline(intercept = 0, slope = 1) + theme_bw()
    idx1 = which(logLRmix > 2 & abs(z) > 2)
	idx2 = which(abs(z_std_diff)>3)
	####Note we made the revision here to include SNPs that are not marginally significant and in high LD with GWAS signals
	idx2_report_index<-which(unlist(lapply(idx2,function(id){
		identifier = ifelse(any(abs(z[which(abs(R[id,])>0.8)])>abs(qnorm(2.5e-8))),1,0)
	}))==1)
	
	idx = union(idx1, idx2[which(idx2_report_index==1)])
	
    if (length(idx) > 0) {
        p = p + geom_point(data = res[idx, ], aes(y = z, x = condmean), 
            col = "red")
    }
    return(list(plot = p, conditional_dist = res))
}

estimate_s_rss<-function (z, R, r_tol = 1e-08, method = "null-mle") 
{
    if (is.null(attr(R, "eigen"))) 
        attr(R, "eigen") = eigen(R, symmetric = TRUE)
    eigenld = attr(R, "eigen")
    if (any(eigenld$values < -r_tol)) 
        warning("The matrix R is not positive semidefinite. Negative ", 
            "eigenvalues are set to zero")
    eigenld$values[eigenld$values < r_tol] = 0
    if (method == "null-mle") {
        negloglikelihood = function(s, z, eigenld) 0.5 * sum(log((1 - 
            s) * eigenld$values + s)) + 0.5 * sum(z * eigenld$vectors %*% 
            ((t(eigenld$vectors) * (1/((1 - s) * eigenld$values + 
                s))) %*% z))
        s = optim(0.5, fn = negloglikelihood, z = z, eigenld = eigenld, 
            method = "Brent", lower = 0, upper = 1)$par
    }
    else if (method == "null-partialmle") {
        colspace = which(eigenld$values > 0)
        if (length(colspace) == length(z)) 
            s = 0
        else {
            znull = crossprod(eigenld$vectors[, -colspace], z)
            s = sum(znull^2)/length(znull)
        }
    }
    else if (method == "null-pseudomle") {
        pseudolikelihood = function(s, z, eigenld) {
            precision = eigenld$vectors %*% (t(eigenld$vectors) * 
                (1/((1 - s) * eigenld$values + s)))
            postmean = rep(0, length(z))
            postvar = rep(0, length(z))
            for (i in 1:length(z)) {
                postmean[i] = -(1/precision[i, i]) * precision[i, 
                  -i] %*% z[-i]
                postvar[i] = 1/precision[i, i]
            }
            return(-sum(dnorm(z, mean = postmean, sd = sqrt(postvar), 
                log = TRUE)))
        }
        s = optim(0.5, fn = pseudolikelihood, z = z, eigenld = eigenld, 
            method = "Brent", lower = 0, upper = 1)$par
    }
    else stop("The method is not implemented")
    return(s)
}

#######################################################################################################################################
#
#
#     Organize GWAS and LD from multiple ancestries to provide summary statistics input for MESuSiE
#
#
#######################################################################################################################################

#' Organize GWAS data from multiple ancestries
#'
#' This function processes and organizes GWAS from two ancestries. It ensures that the dataframes have
#' the necessary columns, checks and warns about any issues, renames columns to a standard format,
#' and outputs a list containing the processed dataframes.
#'
#' @param a1 A GWAS dataframe containing SNP information. At a minimum, the dataframe 
#'           must either have columns "SNP", "beta", "se", and "n" or columns "SNP", "Z", and "n". 
#'           Additional columns like "ref", "alt", "chr", and "pos" are recommended but not mandatory.
#' @param a2 A second GWAS dataframe with SNP details. Similarly to `a1`, it must 
#'           either have columns "SNP", "beta", "se", and "n" or columns "SNP", "Z", and "n". 
#'           Additional columns like "ref", "alt", "chr", and "pos" are recommended but not mandatory.
#' @param element_names A character vector of length 2, naming the resulting list elements. Default is c("data1", "data2").
#'
#' @return A list containing two processed and organized genetic dataframes with user-provided element names.
#'
#' @details
#' The function processes the following:
#' - Checks for required columns.
#' - Computes Z, Beta, and SE values.
#' - Checks and warns for strand-ambiguous and multi-allelic SNPs.
#' - Renames columns to a standard format.
#' - Warns if SNP sets do not match and subsets to common SNPs.
#' - Reports NA values in datasets.
#' - Checks if REF and ALT alleles match and adjusts data if needed.
#' - Suggests providing chromosome and base pair position for further visualization.
#' @export
organize_gwas <- function(a1, a2, element_names=c("data1", "data2")) {
  
  process_data <- function(df) {
    # Normalize column names to lowercase for internal processing
    colnames(df) <- tolower(colnames(df))
    
    # Check if snp column exists
    if (!("snp" %in% colnames(df))) {
      stop("The dataframe does not contain an SNP column.")
    }
    
    # Check required columns
    has_beta_se_N <- all(c("beta", "se", "n") %in% colnames(df))
    has_Z_N <- all(c("z", "n") %in% colnames(df))
    
    if (!has_beta_se_N && !has_Z_N) {
      stop("DataFrame does not contain the required columns")
    }
    
    # Compute Z if beta and se are present
    if (has_beta_se_N && !has_Z_N) {
      df$z <- df$beta / df$se
    }
    
    # Compute se and beta if Z and N are present
    if (has_Z_N && !has_beta_se_N) {
      df$se <- sqrt(1/df$n)
      df$beta <- df$z * df$se
    }
    # Check for strand-ambiguous and multi-allelic SNPs if REF and ALT are provided
	if (all(c("ref", "alt") %in% colnames(df))) {
	  
	  # Identify strand-ambiguous SNPs
	  strand_ambiguous <- df$ref %in% c("A", "T") & df$alt %in% c("A", "T") | 
						  df$ref %in% c("C", "G") & df$alt %in% c("C", "G")
	  if (any(strand_ambiguous)) {
		message("Warning: Detected strand-ambiguous SNPs. Consider reviewing or removing them.")
	  }
	  
	  # Identify multi-allelic SNPs
	  multi_allelic <- nchar(df$ref) > 1 | nchar(df$alt) > 1
	  if (any(multi_allelic)) {
		message("Warning: Detected multi-allelic SNPs. Consider reviewing or removing them.")
	  }
	  
	} else {
	  message("It is suggested to provide REF and ALT allele information, and make sure the REF and ALT alleles are matched up across ancestries. Flip the sign of Beta/Z-score if not matched up.")
	}

    
    # Rename the columns
    colnames(df)[colnames(df) == "snp"] <- "SNP"
    colnames(df)[colnames(df) == "beta"] <- "Beta"
    colnames(df)[colnames(df) == "se"] <- "Se"
    colnames(df)[colnames(df) == "z"] <- "Z"
    colnames(df)[colnames(df) == "n"] <- "N"
	
    # Safely rename REF and ALT if they exist
    if ("ref" %in% colnames(df)) {
      colnames(df)[colnames(df) == "ref"] <- "REF"
    }
    if ("alt" %in% colnames(df)) {
      colnames(df)[colnames(df) == "alt"] <- "ALT"
    }
  
    return(df)
  }
  
  # Check if SNP column exists in both datasets
  if (!("SNP" %in% colnames(a1)) || !("SNP" %in% colnames(a2))) {
    stop("One or both dataframes do not contain an SNP column.")
  }
  
  
	# Process both datasets
	a1 <- process_data(a1)
	a2 <- process_data(a2)
	
	# Record initial SNP counts
	initial_snps_a1 <- nrow(a1)
	initial_snps_a2 <- nrow(a2)

	# Check if SNPs match in both datasets
	common_snps <- intersect(a1$SNP, a2$SNP)

	if (length(common_snps) < nrow(a1) || length(common_snps) < nrow(a2)) {
	warning("The SNPs in the datasets do not match. Subsetting to the common set of SNPs.")
	a1 <- a1[a1$SNP %in% common_snps, ]
	a2 <- a2[a2$SNP %in% common_snps, ]
	}

	# Check and report NAs in datasets
	na_cols_a1 <- colnames(a1)[sapply(a1, function(col) any(is.na(col)))]
	na_cols_a2 <- colnames(a2)[sapply(a2, function(col) any(is.na(col)))]

	if (length(na_cols_a1) > 0) {
	message(paste("Warning: Detected NA values in the first dataset for columns:", paste(na_cols_a1, collapse=", ")))
	message(paste("Counts of NA values for each column in the first dataset:", 
			  paste(paste(na_cols_a1, sapply(a1[na_cols_a1], function(col) sum(is.na(col))), sep=": "), collapse=", ")))
	}

	if (length(na_cols_a2) > 0) {
	message(paste("Warning: Detected NA values in the second dataset for columns:", paste(na_cols_a2, collapse=", ")))
	message(paste("Counts of NA values for each column in the second dataset:", 
			  paste(paste(na_cols_a2, sapply(a2[na_cols_a2], function(col) sum(is.na(col))), sep=": "), collapse=", ")))
	}
	
	# Check if REF and ALT alleles match in both datasets
	if (all(c("REF", "ALT") %in% colnames(a1)) && all(c("REF", "ALT") %in% colnames(a2))) {
	flip_indices <- which(a1$REF == a2$ALT & a1$ALT == a2$REF)
	if(length(flip_indices) > 0) {
	  a2$Beta[flip_indices] <- -a2$Beta[flip_indices]
	  a2$Z[flip_indices] <- -a2$Z[flip_indices]
	  tmp <- a2$REF[flip_indices]
	  a2$REF[flip_indices] <- a2$ALT[flip_indices]
	  a2$ALT[flip_indices] <- tmp
	  warning("The REF and ALT alleles between the datasets are swapped. Adjusting the data by flipping the signs of Beta and Z for consistency.")
	} else if(!identical(a1$REF, a2$REF) || !identical(a1$ALT, a2$ALT)) {
	  stop("The REF and ALT alleles between the datasets are inconsistent and are not simply swapped. Please verify the allele information.")
	}
	}

	# Check if chr and pos are present
	if (!(all(c("chr", "pos") %in% colnames(a1)) && all(c("chr", "pos") %in% colnames(a2)))) {
	message("It is suggested to provide chromosome and base pair position for visualization using MESuSiE_PIP_Plot.")
	}
	# Report
	message(paste("Number of SNPs in the first dataset before organizing:", initial_snps_a1))
	message(paste("Number of SNPs in the second dataset before organizing:", initial_snps_a2))
	message(paste("Number of common SNPs after organizing:", length(common_snps)))
  
	
	# Combine into a list with user-provided element names
	result_list <- setNames(list(a1, a2), element_names)
	
	return(result_list)
}



#' Organize LD matrices and match them to GWAS data.
#'
#' This function checks and organizes two LD matrices to match them with the GWAS data provided by `organize_data`. The LD matrices are adjusted for symmetry if necessary. The SNPs in the column names of the LD matrices should exactly match the SNPs in the GWAS datasets.
#'
#' @param ld1 Matrix: The first LD matrix.
#' @param ld2 Matrix: The second LD matrix.
#' @param gwas_data List: The output from the `organize_data` function containing two GWAS datasets.
#' 
#' @return A list of two LD matrices with names matching the GWAS datasets from `gwas_data`.
#' @export
organize_ld <- function(ld1, ld2, gwas_data) {
  # Ensure gwas_data is the output of organize_data
  if (length(names(gwas_data)) != 2) {
    stop("gwas_data should be the output of organize_data with two elements")
  }
  
  # Function to check and adjust symmetry
  ensure_symmetry <- function(matrix) {
    if (!identical(matrix, t(matrix))) {
      matrix <- (matrix + t(matrix)) / 2
    }
    return(matrix)
  }
  
  # Adjust LD matrices if they aren't symmetric
  ld1 <- ensure_symmetry(ld1)
  ld2 <- ensure_symmetry(ld2)
  
  # Check for NA values in LD matrices
  if (any(is.na(ld1))) {
    stop("The first LD matrix contains NA values. Please rectify before proceeding.")
  }
  if (any(is.na(ld2))) {
    stop("The second LD matrix contains NA values. Please rectify before proceeding.")
  }
  
  # Check if LD matrices have column names
  if (is.null(colnames(ld1)) || is.null(colnames(ld2))) {
    stop("Please provide SNP names for both LD matrices in the column names, and ensure LD matrices are matched up with GWAS data.")
  }
  
  # Check column names of LD matrices with SNP columns of GWAS data
  if (!identical(colnames(ld1), gwas_data[[1]]$SNP)) {
    stop("SNPs in the column names of the first LD matrix do not exactly match or have a different order compared to the SNPs in the first GWAS dataset.")
  }
  if (!identical(colnames(ld2), gwas_data[[2]]$SNP)) {
    stop("SNPs in the column names of the second LD matrix do not exactly match or have a different order compared to the SNPs in the second GWAS dataset.")
  }
  
  # Combine LD matrices into a list with same names as gwas_data
  ld_list <- setNames(list(ld1, ld2), names(gwas_data))
  
  return(ld_list)
}

