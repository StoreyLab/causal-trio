#' auROC
#' 
#' This function calculates the area under the ROC curve.
#'
#'
#' @param pvals A vector of p values.
#' @param causal_indices A vector of causal loci indicies.
#' @param curve A logistic value required by \code{PRROC::roc.curve}.
#' 
#' @importFrom PRROC roc.curve
#'
#' @return A numeric value of auROC.
calc_auroc_pvalues <- function(pvals, causal_indices, curve = FALSE) {
  if ( missing( pvals ) )
    stop( '`pvals` is required!' )
  if ( missing( causal_indices ) )
    stop( '`causal_indices` is required!' )
  if ( length( causal_indices ) == 0 )
    stop( '`causal_indices` must have at least one index!' )
  
  # in some cases there is nothing to do (LMM has singular information matrix)
  # NA is best value to return in that case (scalar)
  if ( is.null( pvals ) )
    return( NA )
  
  # sometimes p-values are missing (GCATest!), treat as worst possible value (p=1)
  if ( anyNA( pvals ) )
    pvals[ is.na( pvals ) ] <- 1
  
  # check range of data here, to complain if it was bad
  if ( any( pvals < 0 ) )
    stop( 'Input p-values included negative values!' )
  if ( any( pvals > 1 ) )
    stop( 'Input p-values included values exceeding 1!' )
  
  # turn p-values into "scores" (for some reason this is needed to have the correct PR curve; is it just the order flip?)
  scores <- - log( pvals )
  
  # separate scores for both classes
  scores_alt <- scores[ causal_indices ]
  scores_nul <- scores[ -causal_indices ]
  
  # make sure both lists are non-empty
  # scores_alt was tested through testing causal_indices
  # just test the other one
  if ( length( scores_nul ) == 0 )
    stop( 'There were no null (non-causal) cases!' )
  
  # generate data, skip curve (default)
  roc <- PRROC::roc.curve( scores_alt, scores_nul, curve = curve )
  
  # return either the full object for plotting, or the AUC only (default)
  return( if (curve) roc else roc$auc)
}


#' auPRC
#' 
#' This function calculates the area under the precision-recall curve.
#'
#'
#' @param pvals A vector of p values.
#' @param causal_indices A vector of causal loci indicies.
#' @param curve A logistic value required by \code{PRROC::pr.curve}.
#' 
#' @importFrom PRROC pr.curve
#'
#' @return A numeric value of auPRC.
calc_auprc_pvalues <- function(pvals, causal_indices, curve = FALSE) {
  if ( missing( pvals ) )
    stop( '`pvals` is required!' )
  if ( missing( causal_indices ) )
    stop( '`causal_indices` is required!' )
  if ( length( causal_indices ) == 0 )
    stop( '`causal_indices` must have at least one index!' )
  
  # in some cases there is nothing to do (LMM has singular information matrix)
  # NA is best value to return in that case (scalar)
  if ( is.null( pvals ) )
    return( NA )
  
  # sometimes p-values are missing (GCATest!), treat as worst possible value (p=1)
  if ( anyNA( pvals ) )
    pvals[ is.na( pvals ) ] <- 1
  
  # check range of data here, to complain if it was bad
  if ( any( pvals < 0 ) )
    stop( 'Input p-values included negative values!' )
  if ( any( pvals > 1 ) )
    stop( 'Input p-values included values exceeding 1!' )
  
  # turn p-values into "scores" (for some reason this is needed to have the correct PR curve; is it just the order flip?)
  scores <- - log( pvals )
  
  # separate scores for both classes
  scores_alt <- scores[ causal_indices ]
  scores_nul <- scores[ -causal_indices ]
  
  # make sure both lists are non-empty
  # scores_alt was tested through testing causal_indices
  # just test the other one
  if ( length( scores_nul ) == 0 )
    stop( 'There were no null (non-causal) cases!' )
  
  # generate data, skip curve (default)
  pr <- PRROC::pr.curve( scores_alt, scores_nul, curve = curve )
  
  # return either the full object for plotting, or the AUC only (default)
  return( if (curve) pr else pr$auc.integral )
}




