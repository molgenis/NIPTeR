#'Calculate chromosomal fraction
#'
#'@param nipt_sample NIPTSample to retrieve chromosomal fraction for
#'@export
chrfractions <- function(nipt_sample)
{
  if(is.null(attr(nipt_sample, "class"))){
    print("Object is not of type Sample")
  }
  else  UseMethod("chrfractions", nipt_sample)
}
#'Get control chromosomes names
#'
#'@param nipt_sample A sample to check wether combined or separated strands 
#'are used
#'@param control_chromosomes Vector with control chromosomes
#'@export
getcontrolchromosomes <- function(nipt_sample, control_chromosomes = control_chromosomes)
{
  if(is.null(attr(nipt_sample, "class"))){
    print("Object is not of type Sample")
  }
  else  UseMethod("getcontrolchromosomes", nipt_sample)
}
#' Get all chromosomal fractions of a control group
#' 
#'@param nipt_control_group The NIPTControlGroup to retrieve the 
#'chromosomal fraction of every autosome for
#'
#'@export
getfractionscontrolgroup <- function(nipt_control_group)
{
  if(is.null(attr(nipt_control_group, "class"))){
    print("Object is not of type NIPT Control Group")
  }
  else  UseMethod("getfractionscontrolgroup", nipt_control_group)
}
#'Get reads per chromosome per control group sample
#'
#'@param nipt_control_group The control group to retrieve reads for
#'@export
getreadscontrolgroup <- function(nipt_control_group)
{
  if(is.null(attr(nipt_control_group, "class"))){
    print("Object is not of type NIPT Control Group")
  }
  else  UseMethod("getreadscontrolgroup", nipt_control_group)
}
#'@export
chrfractions.SeparatedStrands <- function(nipt_sample){
  sapply(X = nipt_sample[[autosomal_chromosome_reads]], FUN = function(x) (rowSums(x) / sum(x)) / 2)
}
#'@export
chrfractions.CombinedStrands <- function(nipt_sample){
  sapply(X = nipt_sample[[autosomal_chromosome_reads]], FUN = function(x) (rowSums(x) / sum(x)))
}

chrreads <- function(nipt_sample){
  sapply(X = nipt_sample[[autosomal_chromosome_reads]], FUN = function(x) rowSums(x)) 
}
#'Retrieve the chromosomal fractions of a chromosome of interest
#'
#'@param nipt_sample NIPTSample to check wether the strands are combined or 
#'separated
#'@param chromo_focus The chromosome of interest
#'@param chromosomal_fracs The chromosomal fractions to extract the 
#'chromosome of interest from
#'@export
retrieve_fractions_of_interest <- function(nipt_sample, chromo_focus, chromosomal_fracs){
  if(is.null(attr(nipt_sample, "class"))){
    print(" ")
  }
  else UseMethod("retrieve_fractions_of_interest", nipt_sample)
}
#'@export
retrieve_fractions_of_interest.CombinedStrands <- function(nipt_sample, chromo_focus, chromosomal_fracs){
  chromosomal_fracs[as.character(chromo_focus),] 
}

#'@export
retrieve_fractions_of_interest.SeparatedStrands <- function(nipt_sample, chromo_focus, chromosomal_fracs){
  chromosomal_fracs[paste0(chromo_focus, "F"),] + chromosomal_fracs[paste0(chromo_focus, "R"), ]
}

setrownamesmatrix <- function(nipt_matrix){
  if (nrow(nipt_matrix) == 22){
    rownames(nipt_matrix) <- rownames_combined_autosomal
  }
  if (nrow(nipt_matrix) == 44){
    rownames(nipt_matrix) <- c(rownames_separated_forward_autosomal, rownames_separated_reverse_autosomal)
  }
  return (nipt_matrix)
}
getsamplenames <- function(nipt_sample){
  nipt_sample[[sample_name]]
}

getcorrectionstatus <- function(nipt_sample, status_type){
  nipt_sample[[status_type]]
}

getstrandtype <- function(nipt_sample){
  class(nipt_sample)[2]
}
#'@export
getcontrolchromosomes.SeparatedStrands <- function(nipt_sample, control_chromosomes = control_chromosomes){
  c(paste0(control_chromosomes, "F"), paste0(control_chromosomes, "R"))
}
#'@export
getcontrolchromosomes.CombinedStrands <- function(nipt_sample, control_chromosomes = control_chromosomes){
  control_chromosomes
}
splitchromosomes <- function(read_counts, chromosomes){
  read_counts[chromosomes,]
}
appendchromosomes <- function(autosomal_reads, sex_reads){
  rbind(autosomal_reads, sex_reads)
}
#'@export
getfractionscontrolgroup.SeparatedStrands <- function(nipt_control_group){
  fraction_table <- sapply(nipt_control_group[[samples]], chrfractions)
  colnames(fraction_table) <- sapply(nipt_control_group[[samples]], getsamplenames)
  rownames(fraction_table) <- c(rownames_separated_forward_autosomal, rownames_separated_reverse_autosomal)
  return(fraction_table)
}
#'@export
getfractionscontrolgroup.CombinedStrands <- function(nipt_control_group){
  fraction_table <- sapply(nipt_control_group[[samples]], chrfractions)
  colnames(fraction_table) <- sapply(nipt_control_group[[samples]], getsamplenames)
  rownames(fraction_table) <- as.character(autosomal_chromosomes)
  return(fraction_table)
}
#'@export
getreadscontrolgroup.SeparatedStrands <- function(nipt_control_group){
  reads_table <- sapply(nipt_control_group[[samples]], chrreads)
  colnames(reads_table) <- sapply(nipt_control_group[[samples]], getsamplenames)
  rownames(reads_table) <- c(rownames_separated_forward_autosomal, rownames_separated_reverse_autosomal)
  return(reads_table)
}
#'@export
getreadscontrolgroup.CombinedStrands <- function(nipt_control_group){
  reads_table <- sapply(nipt_control_group[[samples]], chrreads)
  colnames(reads_table) <- sapply(nipt_control_group[[samples]], getsamplenames)
  rownames(reads_table) <- as.character(autosomal_chromosomes)
  return(reads_table)
}
sumfandrautosomal <- function(nipt_sample){
  summed_reads <- Reduce("+", nipt_sample[[autosomal_chromosome_reads]])
  rownames(summed_reads) <- autosomal_chromosomes
  return (summed_reads)
}
sumfandrsex <- function(nipt_sample){
  summed_reads <- Reduce("+", nipt_sample[[sex_chromosome_reads]])
  rownames(summed_reads) <- XY
  return (summed_reads)
}
testtrainset <- function(nipt_control_group, size_of_train_set){
  indices <- sample(x = 1:length(nipt_control_group[[samples]]), round(size_of_train_set * length(nipt_control_group[[samples]])))
  control_group_train <- as_control_group(nipt_samples = nipt_control_group[[samples]][indices])
  control_group_test <- as_control_group(nipt_samples = nipt_control_group[[samples]][-indices])
  return (list(train_set = control_group_train, test_set = control_group_test))
}
