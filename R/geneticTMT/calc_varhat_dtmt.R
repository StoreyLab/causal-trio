# Transmission-Mean Test Statistic Sampling Variance Estimate
#'
#' This function computes the sampling variance estimate for the TMT statistic 
#' given trio genotypes and children's phenotype. The phenotype can be either 
#' quantitative or dichotomous.
#'
#' @param Y The n-length vector of children' phenotype.
#' @param is_phased Label for unphased/phased genotype. If `FALSE` (default), the
#'     child's genotype is unphased. Else, the child's genotype is phased.
#' @param trio_geno A \eqn{3 \times n}{3 x n} matrix of trio genotype if `FALSE`
#'     for is_phased. Or a \eqn{4 \times n}{4 x n} matrix of trio genotype if
#'     `TRUE` for is_phased. Must contain only values \eqn{(0, 1, 2)}{0, 1, 2}.
#'     The first row is assumed maternal, the second row is assumed paternal, and
#'     the third (or together with the fourth) row is assumed child's genotypes.
#'
#' @return A numeric value of the sampling variance estimate.
#' 
#' @note Undesired dimensions of the genotype matrices will result in errors.
#'
#' @examples 
#' trio_example <- matrix(c(1,1,2,0, 1,2,1,0, 1,1,1,0), nrow=3, byrow=TRUE) 
#' y_quant <- c(11.8,1.7,3.6,12.2)
#' varhat <- calc_varhat_dtmt(y_quant, trio_example)
#' 
#' @export 
calc_varhat_dtmt <- function(Y, trio_geno, is_phased=FALSE){
  if(is_phased){
    geno_c <- trio_geno[3,]+trio_geno[4,]
  } else {
    geno_c <- trio_geno[3,]
  }
  # method 1 for estimating gamma #####################################
  id_het_m <- (trio_geno[1,]==1)
  id_het_p <- (trio_geno[2,]==1)
  id_N <- (id_het_m | id_het_p)
  # calculate var dtmt ########################################################
  # construct assignment vectors
  list_assign <- assignment_index(trio_geno, is_phased)
  W0 <- list_assign$W0
  W1 <- list_assign$W1
  # estimate mu to center trait value
  mu_hat <- calc_mu_hat(Y, W0, W1)
  # sizes of treatment and control
  N <- sum(W1) + sum(W0)
  # calculate the tmt estimand
  dtmt <- ifelse(N>0, 2/N*(sum((Y-mu_hat)*W1) - sum((Y-mu_hat)*W0)), NaN)
  # imbens-rubin approach ######################################################
  T0 <- (W0==1) & (W1==0)
  T1 <- (W0==0) & (W1==1)
  T00 <- (W0==2) & (W1==0)
  T11 <- (W0==0) & (W1==2)
  y <- as.numeric(Y)
  # variance estimates under imbens-rubin approach #############################
  mu_hat_0 <- sum(y[T0]) / sum(T0)
  mu_hat_1 <- sum(y[T1]) / sum(T1)
  mu_hat_00 <- 2*sum(y[T00]) / sum(T00)
  mu_hat_11 <- 2*sum(y[T11]) / sum(T11)
  sigma2_hat_0 <- sum((y[T0]-mu_hat_0)**2) / (sum(T0)-1)
  sigma2_hat_1 <- sum((y[T1]-mu_hat_1)**2) / (sum(T1)-1)
  sigma2_hat_00 <- sum((2*y[T00]-mu_hat_00)**2) / (sum(T00)-1)
  sigma2_hat_11 <- sum((2*y[T11]-mu_hat_11)**2) / (sum(T11)-1)
  varhat <- 4/(N**2)*(sum(T0)*sigma2_hat_0 + sum(T1)*sigma2_hat_1 + sum(T00)*sigma2_hat_00 + sum(T11)*sigma2_hat_11)
  return(varhat)
}
