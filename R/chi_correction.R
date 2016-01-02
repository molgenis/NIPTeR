mean_correction = function(m, bins_overall_mean, reads){
  cor <- bins_overall_mean / sum(m)
  reads * cor
}

sumchiscores = function(bins_list, bins_correct){
  sample_iterations <- 1:length(bins_list)
  bins_overall_mean = sum(as.numeric(unlist(bins_list))) / length(bins_list)
  bins_list_scaled  <- lapply(sample_iterations, (function(x) mean_correction(bins_list[[x]], 
                                                                              bins_overall_mean, bins_correct[[x]])))
  bins_sum_corrected = Reduce("+", bins_list_scaled)
  bins_scaled_expected = bins_sum_corrected / length(bins_list)
  bins_chi_score = bins_list_scaled
  for (i in 1:length(bins_list)) {
    bins_chi_score[[i]] = (bins_scaled_expected - bins_chi_score[[i]])^2 / bins_scaled_expected
  }
  return(Reduce("+", bins_chi_score))
}
correctsamples <- function(nipt_sample, chisumbins, degrees_of_freedom, chi_cutoff, xy_bins){
  construct_sample(autosomal_reads = lapply(X = nipt_sample[[autosomal_chromosome_reads]], FUN = correctbins, 
                                            chisumbins = chisumbins, degrees_of_freedom = degrees_of_freedom, 
                                            chi_cutoff = chi_cutoff), 
                   sex_reads = nipt_sample[[sex_chromosome_reads]], name = nipt_sample[[sample_name]], 
                   correction_status_autosomal = c(nipt_sample[[correction_status_autosomal_chromosomes]], ChiCorrected),
                   correction_status_sex = nipt_sample[[correction_status_sex_chromosomes]])
  
}
correctsamples_include_xy <- function(nipt_sample, chisumbins, degrees_of_freedom, chi_cutoff, chisumbins_xy){
  construct_sample(autosomal_reads = lapply(X = nipt_sample[[autosomal_chromosome_reads]], FUN = correctbins, 
                                            chisumbins = chisumbins, degrees_of_freedom = degrees_of_freedom, 
                                            chi_cutoff = chi_cutoff), 
                   sex_reads = lapply(X = nipt_sample[[sex_chromosome_reads]], FUN = correctbins, chisumbins = chisumbins_xy,
                                      degrees_of_freedom = degrees_of_freedom, chi_cutoff = chi_cutoff),
                   name = nipt_sample[[sample_name]], 
                   correction_status_autosomal = c(nipt_sample[[correction_status_autosomal_chromosomes]], ChiCorrected),
                   correction_status_sex = c(nipt_sample[[correction_status_sex_chromosomes]], ChiCorrected))
}
correctbins = function(bins, chisumbins, degrees_of_freedom, chi_cutoff){
  chi_sum_bins_normalized <- (chisumbins - degrees_of_freedom) / (sqrt( 2 * degrees_of_freedom))

  chi_sum_bins_correction_factor = as.matrix(chisumbins / degrees_of_freedom)
  index = which(chi_cutoff < chi_sum_bins_normalized) 
  bins[index] <- bins[index] / chi_sum_bins_correction_factor[index] 
  return(bins) 
}
#'Performs chi-square based variation reduction
#'
#'@param nipt_sample The NIPTSample object that is the focus of the analysis
#'@param nipt_control_group The NIPTControlGroup object used in the analysis
#'@param chi_cutoff The Z-score cutoff. If a bin has a Z-score above this threshold,
#'it will be corrected
#'@param include_XY Also apply correction to X and Y chromosomes? 
#'
#'@details The chi-squared based variation reduction identifies overdispersed bins within 
#'the control group and corrects these bins in both the sample of interest and the control group. 
#'The function takes in a `NIPTSample` and a `NIPTControlGroup` object, both to be corrected. 
#'For every corresponding bin in the control group a chi-squared score is calculated and this 
#'total score is converted to a normal distribution. Corresponding bins with a normalized score 
#'above _chi_cutoff_ (default 3.5) are corrected by dividing the number of reads by the total 
#'chi-squared score divided by degrees of freedom
#'@return Named list of length 2. The corrected nipt_sample is in index 1 and the 
#'corrected control group in index 2
#'to extract the corrected sample use \code{$sample} or \code{[[1]]}. 
#'To extract the control group from the list use
#'\code{$control_group} or \code{[[2]]}
#'
#'@examples 
#' \dontrun{
#' ##Apply chi-squared based variation reduction method
#' chi_corrected_data <- chicorrect(nipt_sample = gc_LOESS_corrected_sample, 
#'                                  nipt_control_group = subset_loess_corrected_control_group)
#' ##Extract sample and control group
#' loess_chi_corrected_sample <- chi_corrected_data$sample
#' subset_loess_chi_corrected_control_group <- chi_corrected_data$control_group
#' }
#' 
#'@export
chi_correct <- function(nipt_sample, nipt_control_group, chi_cutoff = 3.5, include_XY = F){
  degrees_of_freedom  = length(nipt_control_group[[samples]]) - 1
  controlbins <- lapply(X = nipt_control_group[[samples]], FUN = function(x) Reduce("+", x[[autosomal_chromosome_reads]]))
  xy_bins <- lapply(X = nipt_control_group[[samples]], FUN = function(x) x[[sex_chromosome_reads]])
  chisumbins <- sumchiscores(bins_list = controlbins, bins_correct = controlbins)
  if (include_XY == T){
    xy_controlbins <- lapply(xy_bins, function(x) Reduce("+", x))
    chisumbins_xy = sumchiscores(bins_list = controlbins, bins_correct = xy_controlbins)
    correctedcontrolsamples <- lapply(X = nipt_control_group[[samples]], FUN = correctsamples_include_xy, chisumbins = chisumbins, 
                                      degrees_of_freedom = degrees_of_freedom, chi_cutoff = chi_cutoff, 
                                      chisumbins_xy = chisumbins_xy)
    correctedsample <- correctsamples_include_xy(nipt_sample = nipt_sample, chisumbins = chisumbins, 
                                                 degrees_of_freedom = degrees_of_freedom, chi_cutoff = chi_cutoff,
                                                chisumbins_xy = chisumbins_xy)
  }
  else{
    correctedcontrolsamples <- lapply(X = nipt_control_group[[samples]], FUN = correctsamples, chisumbins = chisumbins, 
                                      degrees_of_freedom = degrees_of_freedom,
                                      chi_cutoff = chi_cutoff)
    correctedsample <- correctsamples(nipt_sample = nipt_sample, chisumbins = chisumbins, degrees_of_freedom = degrees_of_freedom,
                                      chi_cutoff = chi_cutoff)
  }
  return (list(sample = correctedsample, control_group = as_control_group(correctedcontrolsamples, nipt_control_group[[description]])))
}