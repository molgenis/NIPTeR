
gc_correct_NIPTSample_bin <- function(nipt_object, include_XY, ref_genome ){
  correct_sample_bin(nipt_object, include_XY = include_XY, ref_genome)
}
gc_correct_NIPTControlGroup_bin <- function(nipt_object, include_XY, ref_genome){
  as_control_group(nipt_samples = lapply(nipt_object[[samples]], correct_sample_bin, 
                                         include_XY = include_XY, ref_genome))
}

GetMeanNumberOfReadsGcStep <- function(indices.gc.percentages, sample){
  avg.reads.gc.set <- NULL  
  for (i in 1:length(indices.gc.percentages))  {
    reads <- sample[indices.gc.percentages[[i]]]
    #If total number of reads in GC interval is not 0, determine average number of reads per bin, ignoring bins with no reads
    if(sum(reads != 0))    {
      avg.reads.gc.set[i] <- sum(reads) /length(reads[reads > 0])
    }
  }
  return(avg.reads.gc.set)
}
#Corrects the number of reads in a given GC interval using the weights
CorrectGcWithWeights <- function(sample, indices.gc.percentages, weights){
  for (i in 1:length(indices.gc.percentages)){
    #For a given GC interval, subset matching bins and correct with weight
    sample[indices.gc.percentages[[i]]] <- sample[indices.gc.percentages[[i]]] * weights[i]
  }
  return(sample)
}

correct_sample_bin <- function(nipt_sample, include_XY, ref_genome){
  sample <- Reduce("+", nipt_sample[[autosomal_chromosome_reads]])

  rownames(sample) <- autosomal_chromosomes
  if (ref_genome == "hg37"){
  gc_percentages <- gc_percentages_hg37
  indices.gc.percentages <- indices.gc.percentages.h37
  sex_chromosome_indices <- sex_chromosome_indices_37
  }
  if (ref_genome == "hg38"){
  gc_percentages <- gc_percentages_hg38
  indices.gc.percentages <- indices.gc.percentages.h38
  sex_chromosome_indices <- sex_chromosome_indices_38
  }
  gc_percentages_autosomal <- gc_percentages[autosomal_chromosomes, ]
  #remove reads which have no GC % count
  no_gc_count <- which(gc_percentages_autosomal < 0)
  sample[no_gc_count] <- 0
  #Gets mean number of reads per bin for a GC interval
  avg.reads.gc.interval <- GetMeanNumberOfReadsGcStep(indices.gc.percentages, sample)
  #Calculates mean global average of n. of reads per bin
  n.of.reads <- sample[(unlist(indices.gc.percentages))]
  n.of.bins <- length(n.of.reads[n.of.reads > 0]) 
  global.mean.n.of.bins <- sum(sample) / n.of.bins
    #Calculate weights
  weights <- global.mean.n.of.bins / avg.reads.gc.interval
  weights[which(is.na(weights))] <- 0
  #Correct read counts using weights, forward and reverse apart
  corrected_autosomal <- lapply(nipt_sample[[autosomal_chromosome_reads]], CorrectGcWithWeights, indices.gc.percentages = indices.gc.percentages,
                            weights = weights)
  if (include_XY == TRUE){
    corrected_sex <- lapply(nipt_sample[[sex_chromosome_reads]], CorrectGcWithWeights, 
                                  indices.gc.percentages = sex_chromosome_indices,
                                  weights = weights)
  corrected_sample <- construct_sample(autosomal_reads = corrected_autosomal, sex_reads = corrected_sex, 
                                         name = nipt_sample[[sample_name]], correction_status_autosomal = c(nipt_sample[[correction_status_autosomal_chromosomes]],
                                                                                                            GCcorrected),
                                         correction_status_sex = c(nipt_sample[[correction_status_sex_chromosomes]], GCcorrected))
  }
  else{
  corrected_sample <- construct_sample(autosomal_reads = corrected_autosomal, sex_reads = nipt_sample[[sex_chromosome_reads]], 
                                        name = nipt_sample[[sample_name]], correction_status_autosomal = c(nipt_sample[[correction_status_autosomal_chromosomes]], 
                                                                                                           GCcorrected),
                                        correction_status_sex = c(nipt_sample[[correction_status_sex_chromosomes]]))
  }
  return (corrected_sample)
}