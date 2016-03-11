RetrieveSubsets <- function(candidates, number_of_elements){
  return(set_combn(candidates, m=number_of_elements))
}
CalculateVariation <- function(denominators, chromosomal_frac_control_reads, chr_focus){
  possible_denominators <- unlist(denominators)
  if (length(possible_denominators) == 1){
    mean_subset <- mean(chromosomal_frac_control_reads[chr_focus,] / chromosomal_frac_control_reads[possible_denominators,])
    sd_subset <- stats::sd((chromosomal_frac_control_reads[chr_focus,] / (chromosomal_frac_control_reads[possible_denominators,])))
  }
  else{
    mean_subset <- mean(chromosomal_frac_control_reads[chr_focus,] / (colSums(chromosomal_frac_control_reads[possible_denominators,])))
    sd_subset <- stats::sd(chromosomal_frac_control_reads[chr_focus,] / (colSums(chromosomal_frac_control_reads[possible_denominators,])))
  }
  return(list(sd_subset , mean_subset))
}
#'Prepare NCV calculation
#'
#'
#'Determine the best NCV chromosomes,  calculate NCV scores 
#'and asses normal distribution control group using Shapiro-Wilk test
#'@param nipt_control_group The NIPTControlGroup object used in the analysis
#'@param chr_focus Integer.The chromosome of interest. Most commonly chromosome 13, 18 or 21. 
#'However, every autosomal chromosome can be predicted
#'@param max_elements Integer, The maximum number of denominator chromosomes. 
#'@param exclude_chromosomes Integer. Exclude which autosomal chromosomes as potential predictors? 
#'Default potential trisomic chromosomes 13, 18 and 21 are exluded. 
#'@param include_chromosomes Integer. Which potential trisomic chromosomes
#' (13,18 and 21) to include?
#'@param use_test_train_set Boolean. Use a test and train set?
#'@param size_of_train_set Double The size of the train set expressed in a decimal. 
#'Default is 0.6 (60\% of the control group samples)
#'
#'@details
#The Normalized Chromosome Value or NCV (Sehnert et al., 2011) method selects a subset of 
#'chromosomes to calculate the chromosomal fractions. The 'best' subset is the set which 
#'yields the lowest coefficient of variation for the chromosomal fractions of the chromosome 
#'of interest in the control group. Because a brute force approach is used to determine the 
#'best subset, which can be computationally intensive,this method is divided into two functions, 
#'prepare_ncv and calculate_ncv. prepare_ncv returns a template object (NCVTemplate) for 
#'a given chromosome of interest and the control group used. This template can be used for 
#'any number of analyses. If the control group or chromosome of interest changes, 
#'a new template must be made.
#'
#'
#'The ncv_template object is a list containing:
#'\itemize{
#'\item Character \strong{denominators} The set of denominator chromosomes
#'\item Character \strong{focus_chromosome}The chromosome of interest used for this `NCVTemplate` object
#'\item Character \strong{nipt_sample_names} The sample names of the test set samples
#'\item Character \strong{correction_status} The correction status(es) of the
#'control group samples
#'\item Data.frame \strong{control_group_Z_scores} The NCV scores for the test set samples
#'\item Character \strong{potential_denominators} The total pool of denominators the best denominators
#' are selected from
#'\item Numeric \strong{control_group_statistics} Named num of length 3, the first field being the mean 
#'(name mean), the second field is the standard deviation (name SD) and the third field is the P value 
#'of the Shapiro-Wilk test (name Shapiro_P_value)
#' }
#'If a Test and Train set is used the ncv_template object also includes:
#'\itemize{
#'\item Character \strong{sample_names_train_set} The sample name where the model
#'is trained on
#'\item Numeric \strong{train_set_statistics} Mean, SD and Shapiro-Wilk test
#' P value of the Z scores of the train set
#' \item Data.frame \strong{train_set_Zscores} The Z scores of the train set
#'}
#'
#'
#'
#'@return ncv template object
#'
#'@examples 
#' \dontrun{
#' ##Create NCVTemplates for chromosome 13 with max 9 denominators and default settings, so:
#' ##All autosomals chromosomes are potential predictors, 
#' ##except the potential trisomic chromosomes 13, 18 and 21
#' new_ncv_template_13 <- prepare_ncv(nipt_control_group = control_group, 
#'                                    chr_focus = 13, max_elements = 9)
#' }
#'
#'@references
#'\href{http://www.ncbi.nlm.nih.gov/pubmed/21519036}{Sehnert et al.} 
#'
#'@export
prepare_ncv <- function(nipt_control_group, chr_focus, max_elements, exclude_chromosomes = NULL, include_chromosomes = NULL,
                            use_test_train_set =T, size_of_train_set = 0.6){
  exclude_chromosomes <- as.character(exclude_chromosomes)
  include_chromosomes <- as.character(include_chromosomes)
  min_variation_subset <- list()
  min_variation_vc <- NULL
  chromosomal_frac_control_reads <- sapply(lapply(nipt_control_group[[samples]], sumfandrautosomal), rowSums)
  chromosomal_frac_control_reads_train <- chromosomal_frac_control_reads
  if (use_test_train_set == T){
    indices <- sample(x = 1:length(nipt_control_group[[samples]]), round(size_of_train_set * length(nipt_control_group[[samples]])))
    nipt_control_group_train <- as_control_group(nipt_samples = nipt_control_group[[samples]][indices])
    nipt_control_group <- as_control_group(nipt_samples = nipt_control_group[[samples]][-indices])
    chromosomal_frac_control_reads_train <- sapply(lapply(nipt_control_group_train[[samples]], sumfandrautosomal), rowSums)
    rownames(chromosomal_frac_control_reads) <- as.character(autosomal_chromosomes)
    chromosomal_frac_control_reads <- sapply(lapply(nipt_control_group[[samples]], sumfandrautosomal), rowSums)
  }
  rownames(chromosomal_frac_control_reads_train) <- as.character(autosomal_chromosomes)
  rownames(chromosomal_frac_control_reads) <- as.character(autosomal_chromosomes)
  control_chromosomes <- control_chromosomes[!control_chromosomes %in% c(exclude_chromosomes, chr_focus)]
  
  if (!is.null(include_chromosomes)){
    control_chromosomes <- c(control_chromosomes, include_chromosomes)
  }
  for (n_of_elements in 1:max_elements)  {
    denominators <- as.list(RetrieveSubsets(candidates = control_chromosomes, number_of_elements = n_of_elements))
    variation_list <- lapply(FUN = CalculateVariation, X=denominators,
                             chromosomal_frac_control_reads= chromosomal_frac_control_reads_train, chr_focus=chr_focus)
    variation <- sapply(variation_list, function(x) Reduce("/", x))
    min_variation_subset[[n_of_elements]] <- as.numeric(denominators[[which.min(variation)]])
    min_variation_vc[n_of_elements] <- variation[which.min(variation)]
  }
  denominators <- min_variation_subset[[which.min(min_variation_vc)]]
  ncv_reads <- get_ncv_reads(chromosomal_frac_control_reads = chromosomal_frac_control_reads, chr_focus = chr_focus, 
                             denominators = denominators )
  scores <- ZScoresControl(denominators = denominators, ncv_reads = ncv_reads,
                           chr_focus = chr_focus, samplenames = sapply(nipt_control_group[[samples]], getsamplenames))
  if (use_test_train_set == F){
    new_ncv_template <- ncv_template(denominators = denominators, chromo_focus = chr_focus, 
                                     nipt_sample_names = sapply(nipt_control_group[[samples]], getsamplenames),
                                     correction_status = unique(nipt_control_group[[correction_status_autosomal_chromosomes]]), 
                                     scores = scores, potential_denominators = control_chromosomes,
                                     statistics = c(mean = mean(ncv_reads), SD = stats::sd(ncv_reads), 
                                                    Shapiro_Wilk_P_value = shapiro.test(scores$NCVscore)$p.value),
                                     type = class(nipt_control_group)[2])
  }
  else{
    ncv_reads_train <- get_ncv_reads(chromosomal_frac_control_reads = chromosomal_frac_control_reads_train, chr_focus = chr_focus, 
                                     denominators = denominators )
    scores_train <- ZScoresControl(denominators = denominators, ncv_reads = ncv_reads_train,
                                   chr_focus = chr_focus, samplenames = sapply(nipt_control_group_train[[samples]], getsamplenames))
    new_ncv_template <- ncv_template(denominators = denominators, chromo_focus = chr_focus, 
                                     nipt_sample_names = sapply(nipt_control_group[[samples]], getsamplenames),
                                     correction_status = unique(nipt_control_group[[correction_status_autosomal_chromosomes]]), 
                                     scores = scores, potential_denominators = control_chromosomes,
                                     statistics = c(mean = mean(ncv_reads), SD = stats::sd(ncv_reads), 
                                                    Shapiro_Wilk_P_value = shapiro.test(scores$NCVscore)$p.value),
                                     type = class(nipt_control_group)[2], 
                                     sample_names_train_set = sapply(nipt_control_group_train[[samples]], getsamplenames),
                                     train_set_statistics = c(mean = mean(ncv_reads_train), SD = stats::sd(ncv_reads_train), 
                                                              Shapiro_Wilk_P_value = shapiro.test(scores_train$NCVscore)$p.value),
                                     train_set_Zscores = scores_train)  
  }
  
  
  return(new_ncv_template)
}
get_ncv_reads <- function(chromosomal_frac_control_reads, chr_focus, denominators){
  if (length(denominators) == 1){
    return (ncv_reads <- chromosomal_frac_control_reads[chr_focus, ] / chromosomal_frac_control_reads[denominators,])
  }
  else{
    return (ncv_reads <- chromosomal_frac_control_reads[chr_focus, ] / colSums(chromosomal_frac_control_reads[denominators,]))
  }
}
ZScoresControl <- function(denominators, ncv_reads, chr_focus, samplenames){
  scores <- as.data.frame((ncv_reads - mean(ncv_reads)) / stats::sd(ncv_reads), row.names = samplenames)  
  colnames(scores) <- "NCVscore"
  return(scores)
}
get_ncv_statistics <- function(ncv_reads){
  return(c(mean(ncv_reads), stats::sd(ncv_reads)))
}
