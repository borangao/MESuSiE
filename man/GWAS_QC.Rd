% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/meSuSie_QC.R
\name{GWAS_QC}
\alias{GWAS_QC}
\title{GWAS Quality Control Function}
\usage{
GWAS_QC(dat, MAF_threshold)
}
\arguments{
\item{dat}{A data frame containing the GWAS data. Columns "CHR", "POS", "CHR_POS", 
"REF", "ALT", and (optionally) "MAF" are required.}

\item{MAF_threshold}{The threshold for Minor Allele Frequency filtering. SNPs 
with MAF below this threshold or above (1 - threshold) will be removed.}
}
\value{
A filtered data frame after performing the quality control steps.
}
\description{
This function performs several quality control steps on GWAS data:
1. Removes strand ambiguous variants.
2. Removes multi-allelic variants.
3. Excludes the MHC complex region.
4. Filters SNPs based on Minor Allele Frequency (MAF).
}
