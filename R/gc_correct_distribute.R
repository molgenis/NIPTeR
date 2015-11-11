GetCorrectionFactor <- function(sample, span, gc_percentages, median_autosomals){
  fit.loess <-loess(sample ~ gc_percentages, span = span)
  fitted_values <- cbind(fit.loess$x, median_autosomals /  fit.loess$fitted)
  return(list(median_autosomals / fit.loess$fitted, fitted_values))
}
correct.reads <- function(sample, correction.factor, indices){
  sample[indices] <- sample[indices] * correction.factor
  return(sample)
}
correct_reads_sex_chromosome <- function(bins, gc_percentages,fitted_values){
  for (chromosome in 1:nrow(bins))
  {
    gc_bins <- unname(which(gc_percentages[22+chromosome,] > 0))
    for (bin in gc_bins){
      bins[chromosome, bin] <- bins[chromosome, bin] * (fitted_values[order(abs(gc_percentages[(22+chromosome),bin] - fitted_values[,1]))[1],][2])
    }
  }
  return(bins)
}

gc_correct_NIPTSample_loess <- function(nipt_object, span = 0.75, include_XY = include_XY, ref_genome){
  corrected_sample <- correct_sample_loess(nipt_object, span = 0.75, include_XY = include_XY, 
                                           ref_genome = ref_genome)
}
gc_correct_NIPTControlGroup_loess <- function(nipt_object, span = 0.75, include_XY = include_XY, ref_genome){
  as_control_group(nipt_samples = lapply(nipt_object[[samples]], correct_sample_loess, include_XY = include_XY,
                                         ref_genome = ref_genome))
}

correct_sample_loess <- function(nipt_sample, span = 0.75, include_XY = include_XY, ref_genome){
 
  sample <- Reduce("+", nipt_sample[[autosomal_chromosome_reads]])
  rownames(sample) <- autosomal_chromosomes
  if (ref_genome == "hg37"){
    gc_percentages <- gc_percentages_hg37
  }
  if (ref_genome == "hg38"){
    gc_percentages <- gc_percentages_hg38
  }
  gc_percentages_autosomal <- gc_percentages[autosomal_chromosomes, ]
  indices_autosomals <- (which(gc_percentages_autosomal > 0 & sample))
  median_autosomals <- median(sample[indices_autosomals])
  correction_factor_values <- GetCorrectionFactor(sample = as.vector(sample[indices_autosomals]), span = span, 
                                                  gc_percentages = as.vector(gc_percentages_autosomal[indices_autosomals]), 
                                                  median_autosomals = median_autosomals)
  correction.factor <- correction_factor_values[[1]]
  fitted_values <- correction_factor_values[[2]]
  corrected_autosomal <- lapply(nipt_sample[[autosomal_chromosome_reads]],  correct.reads, correction.factor = correction.factor, 
                                indices = indices_autosomals)
  if (include_XY == T){
    corrected_sex <- lapply(X = nipt_sample[[sex_chromosome_reads]], FUN = correct_reads_sex_chromosome, 
                            gc_percentages = gc_percentages, fitted_values = fitted_values)
    corrected.sample <- construct_sample(autosomal_reads = corrected_autosomal, sex_reads = corrected_sex, 
                                         name = nipt_sample[[sample_name]], correction_status_autosomal = c(nipt_sample[[correction_status_autosomal_chromosomes]], GCcorrected),
                                         correction_status_sex = c(nipt_sample[[correction_status_sex_chromosomes]], GCcorrected))
  }
  else{
  corrected.sample <- construct_sample(autosomal_reads = corrected_autosomal, sex_reads = nipt_sample[[sex_chromosome_reads]], 
                                       name = nipt_sample[[sample_name]], correction_status_autosomal = c(nipt_sample[[correction_status_autosomal_chromosomes]], GCcorrected),
                                       correction_status_sex = c(nipt_sample[[correction_status_sex_chromosomes]]))
  }
  return(corrected.sample)
}
