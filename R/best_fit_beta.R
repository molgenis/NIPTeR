#'Best matching control group by least sum of squares
#' 
#'The matchcontrolgroup function determines how well an NIPTSample fits within the NIPTControlGroup 
#'@param nipt_sample The NIPTSample object that is the focus of the analysis
#'@param nipt_control_group The NIPTControlGroup object used in the analysis
#'@param mode The function mode. This can either be \emph{"subset"} or \emph{"report"}.
#'Mode \emph{"subset"} means the return value will be a new `NIPTControlGroup` object containing  \emph{n}
#'samples. When mode \emph{"report"} is used the output is a matrix containing the sum of squares score 
#'of the differences between the chromosomal fractions of the sample and the control for every control 
#'sample, sorted in increasing score.
#'@param n_of_samples The length of the resulting NIPTControlGroup. Only applicable if mode \emph{"subset"} is used.
#'  
#'@param include_chromosomes integer. Include potential trisomic chromosomes into the comparison? Default = NULL,
#' meaning chromosomes 13, 18 and 21 are not included
#'@param exclude_chromosomes integer.Exclude other autosomal chromosomes besides chromosomes 
#'13, 18 and 21? Default = NULL
#' 
#'@details 
#'The `matchcontrolgroup` function determines how well an NIPTSample fits within the NIPTControlGroup
#'and, if needed, makes a subset `NIPTControlGroup` of length \emph{n}.
#' 
#'@return The output for mode \emph{subset} is a new `NIPTControlGroup` composed of _n_ samples. 
#'The output for mode \emph{report} is a matrix with a single column containing the sum of squares
#'in ascending order.
#' 
#' @examples 
#' \dontrun{
#' ##Mode report
#' scores_control_group <- matchcontrolgroup(nipt_sample = sample_of_interest, 
#'                                           nipt_control_group = control_group, 
#'                                           mode = "report", include_chromosomes = c(13,18))
#'
#' ##Mode subset
#' subset_control_group <- matchcontrolgroup(nipt_sample = sample_of_interest, 
#'                                           nipt_control_group = control_group, 
#'                                           mode = "subset", n_of_samples = 50)
#' } 
#' 
#' 
#' @export
match_control_group <-function(nipt_sample, nipt_control_group, mode, n_of_samples, 
                             include_chromosomes = NULL, exclude_chromosomes = NULL)
{
  if (!is.null(include_chromosomes)){
  control_chromosomes <- control_chromosomes[!control_chromosomes %in% exclude_chromosomes]
  }
  if (!is.null(include_chromosomes)){
    control_chromosomes <- c(control_chromosomes, include_chromosomes)
  }
  nipt_control_group_samples <- nipt_control_group[[samples]]
  control_group_fractions <- cbind(NULL, sapply(X = nipt_control_group_samples , FUN = chrfractions))
  sample_fractions <- sapply(list(nipt_sample), chrfractions)
  fractions_squared <- apply(X = control_group_fractions, MARGIN = 2, FUN = function(x, y) {(x - y)^2}, y = sample_fractions)
  control_chromosomes <- getcontrolchromosomes(nipt_sample, control_chromosomes)
  fractions_squared <- setrownamesmatrix(fractions_squared)
  sum_of_squares <- apply(fractions_squared[control_chromosomes,], 2, sum)
  if (mode == "subset"){
  sorted_scores <- order(sum_of_squares)
  return(as_control_group(nipt_samples = nipt_control_group[[samples]][sorted_scores[1:n_of_samples]], 
                   control_group_type = paste("Fitted to", nipt_sample[[sample_name]])))
  }
  if (mode == "report"){
    names(sum_of_squares) <- sapply(nipt_control_group[[samples]], getsamplenames)
    sorted <- as.matrix(sort(sum_of_squares))
    colnames(sorted) <- "Sum_of_squares"
    return(sorted)
  }
  else{
    stop("Invalid mode selected", call. = F)
  }
}

