#'Calculate 'standard' Z-score
#'
#'@details In the Z-score approach, introduced by Chiu et al in 2008, 
#'the chromosomal fraction of interest of a sample is compared 
#'to the chromosomal fractions of interest of the reference samples, 
#'the `NIPTControlGroup` object.
#'The output of the function is an object of class `ZscoreResult`. It is a named list containing seven fields:
#'\itemize{
#'\item numeric \strong{sample_Zscore} The Z score for the sample of interest for the sample of interest
#'\item named num \strong{control_group_statistics} Named num of length 3, the first field being the 
#'mean (name mean), the second field is the standard deviation (name SD) 
#'and the third field is the P value of the Shapiro-Wilk test (name Shapiro_P_value)
#'\item matrix \strong{control_group_Zscores} containing the Z scores of the chromosome of interest for 
#'all used control samples
#'\item integer \strong{focus_chromosome} The chromosome of interest. Most commonly chromosome 13, 18 or 21.
#' However, every autosomal chromosome can be predicted
#'\item string \strong{control_group_sample_names} The sample names of all control group samples used in the
#' analysis
#'\item string \strong{correction_status} The correction status of the control group
#'\item string \strong{sample_name} The sample_name of the sample of interest
#'}
#'@param nipt_sample  The NIPTSample object that is the focus of the analysis
#'@param nipt_control_group The NIPTControlGroup object used in the analysis
#'@param chromo_focus The chromosome of interest. Most commonly chromosome 13, 18 or 21. 
#'However, every autosomal chromosome can be predicted
#'
#'@return ZscoreResult object
#'
#'@examples 
#' \dontrun{
#' z_score_result_13 <- calculate_z_score(nipt_sample = sample_of_interest, 
#'                                        nipt_control_group = control_group, 
#'                                        chromo_focus = 13)
#' }
#'
#'@export
calculate_z_score <- function(nipt_sample, nipt_control_group, chromo_focus){
  control_group_fractions <- sapply(X = nipt_control_group[[samples]], FUN = chrfractions)
  sample_fractions <- sapply(list(nipt_sample), chrfractions)
  control_group_fractions <- setrownamesmatrix(control_group_fractions)
  sample_fractions <- setrownamesmatrix((sample_fractions))
  chromo_focus_fractions <- getstats(fractions = control_group_fractions, chromo_focus = chromo_focus)
  mean_control <- mean(chromo_focus_fractions)
  sd_control <-  stats::sd(chromo_focus_fractions)
  sample_focus_fractions <- getstats(fractions = sample_fractions, chromo_focus = chromo_focus)
  sample_score <- unname((sample_focus_fractions - mean_control) /sd_control)
  control_scores <- scale(chromo_focus_fractions)
  control_names <- sapply(nipt_control_group[[samples]], getsamplenames)
  dimnames(control_scores) <- list(control_names, "Z_score")
  new_z_result <- z_score_template(statistics = c(mean = mean_control, SD = sd_control, 
                                                  Shapiro_P_value = shapiro.test(control_scores)$p.value),
                                   control_group_scores = control_scores, chromo_focus = chromo_focus,
                                   nipt_sample_names = control_names, correction_status = unique(nipt_control_group$correction_status),
                                   sample_name = getsamplenames(nipt_sample), sample_Zscore = sample_score,
                                   type = class(nipt_control_group)[2])
  return(new_z_result)
}

z_score_template <- function(statistics, control_group_scores, chromo_focus, nipt_sample_names, correction_status,
                             sample_Zscore, sample_name, type){
  new_z_score_template <- list(sample_Zscore = sample_Zscore, control_group_statistics = statistics, 
                               control_group_Zscores= control_group_scores,
                               focus_chromosome = as.character(chromo_focus),
                               control_group_sample_names = nipt_sample_names, 
                               correction_status = correction_status,
                               sample_name = sample_name)
  
  class(new_z_score_template) <- c(Z_template_class, type)
  
  return(new_z_score_template)
  
}

getstats <- function(fractions, chromo_focus){
 
  if (nrow(fractions) == 44){
 
  chromo_focus_fractions <- fractions[chromo_focus,] + fractions[chromo_focus + 22,]
  }
  else{
  chromo_focus_fractions <- fractions[chromo_focus,]
  }
  return(chromo_focus_fractions)
}