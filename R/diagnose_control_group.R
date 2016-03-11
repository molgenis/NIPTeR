#'Diagnose control group
#'
#'Compute a regular Z-score for every chromosome of every sample in a NIPTControlGroup object
#'
#'@param nipt_control_group The NIPTControlGroup object to diagnose
#'
#'@details This function computes a regular Z-score for every chromosome of
#'every sample in a NIPTControlGroup object. It returns a named list with 
#'diagnostics information.
#'
#'The function returns a named list with 3 fields:
#'\itemize{
#'\item \strong{Z_scores} A matrix containing Z-scores for every sample and every chromosome
#'\item \strong{abberant_scores} Dataframe with samplename and chromosome of Z-scores outside -3  3 range 
#'\item \strong{control_group_statistics} Matrix with mean, standard deviation and P value of Shapiro-Wilk test
#'}
#'@return named list
#'
#'@examples 
#' \dontrun{
#' diagnose_control_group(nipt_control_group = control_group)
#' }
#' 
#'@export
diagnose_control_group <- function(nipt_control_group){
  fracs <- sapply(nipt_control_group[[samples]], chrfractions)
  control_group_scores <- t(apply(X = fracs, FUN = scale, MARGIN = 1))
  colnames(control_group_scores) <- sapply(nipt_control_group[[samples]], getsamplenames)
  control_group_scores <- setrownamesmatrix(control_group_scores)
  a <- rbind(which(control_group_scores > 3, arr.ind = T), which(control_group_scores < -3, arr.ind = T))
  abberants <- data.frame(t(apply(X = a, FUN = get_abberant_scores, score_matrix=control_group_scores, MARGIN = 1 )))
  if (length(abberants) > 0){
    rownames(abberants) <- NULL
    colnames(abberants) <- c("Chromosome", "Sample_name", "Z_score")
  }
  else{
    abberants = NULL
  }
  statistics <- t(apply(control_group_scores, MARGIN = 1, FUN = get_mean_sd_shapiro))
  statistics <- setrownamesmatrix(statistics)
  
  return(list(Z_scores = control_group_scores, abberant_scores = abberants, 
              control_group_statistics = statistics))
  
}
get_abberant_scores <- function(row, score_matrix){
  c(dimnames(score_matrix)[[1]][row[1]], dimnames(score_matrix)[[2]][row[2]], score_matrix[row[1], row[2]])
}
get_mean_sd_shapiro <- function(row){
  c(mean = mean(row), SD = stats::sd(row), Shapiro_P_value = shapiro.test(row)$p.value)
}
