#'Use an NCV template to calculate a NCV score for sample of interest
#'
#'@param nipt_sample nipt_sample object of interest
#'@param ncv_template ncv_template object, result from \code{\link{prepare_ncv}}
#'
#'@details
#'\code{\link{prepare_ncv}}
#'
#'
#'@return ncv_result object
#'
#'@examples 
#' \dontrun{
#' ##Use NCVTemplate to get NCV scores for the sample of interest
#' ncv_score_13 <- calculate_ncv_score(nipt_sample = sample_of_interest, 
#'                                     ncv_template = new_ncv_template_13)
#' }
#'
#'@references
#'\href{http://www.ncbi.nlm.nih.gov/pubmed/21519036}{Sehnert et al.} 
#'
#'@export
calculate_ncv_score <- function(nipt_sample, ncv_template){
  reads <- (rowSums(sumfandrautosomal(nipt_sample)))
  normalized_chromosome <- reads[ncv_template$focus_chromosome] / sum(reads[ncv_template$denominators])
  ncv_score <- (normalized_chromosome - ncv_template$control_group_statistics[1]) / ncv_template$control_group_statistics[2]
  
  ncv_template$sample_score <- ncv_score
  ncv_template$sample_name <- nipt_sample$name  
  
  class(ncv_template)[1] <- NCV_result_class
  
  return(ncv_template)
}