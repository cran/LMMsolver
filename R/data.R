#' Simulated QTL mapping data set
#'
#' @format A data.frame with 180 rows and 6 columns.
#' \describe{
#' \item{cross}{Cross ID, two populations, AxB and AxC}
#' \item{ind}{Genotype ID}
#' \item{pA}{Probability that individual has alleles from parent A}
#' \item{pB}{Probability that individual has alleles from parent B}
#' \item{pC}{Probability that individual has alleles from parent C}
#' \item{pheno}{Simulated phenotypic value}
#' }
"multipop"

#' Simulated Biomass as function of time using APSIM wheat.
#'
#' @format A data.frame with 121 rows and 4 columns.
#' \describe{
#' \item{env}{Environment, Emerald in 1993}
#' \item{geno}{Simulated genotype g001}
#' \item{das}{Days after sowing}
#' \item{biomass}{Simulated biomass using APSIM; medium measurement error added}
#' }
#' @references Bustos-Korts et al. (2019) Combining Crop Growth Modeling and
#' Statistical Genetic Modeling to Evaluate Phenotyping Strategies
#' \doi{10.3389/FPLS.2019.01491}
"APSIMdat"
