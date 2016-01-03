#'Convert list of nipt samples to nipt control group
#'
#'
#'@param nipt_samples List of nipt_sample objects to be combined to a control
#'group
#'@param control_group_type Control group type, either 'generic control group' 
#'or 'fitted to sample'. Leave this argument blank
#'
#'@details
#'This function returns an S3 object of class nipt_control_group. It is a list 
#'with 3 items:
#'
#'\itemize{
#'\item List \strong{Samples} nipt_sample objects in the control group
#'\item Character \strong{Correction_status} Correction_status(es) in the 
#'control group
#'\item Character \strong{Samplenames} The sample names of samples present 
#'in the control group
#'}
#'Read count strategy should be uniform in all samples in a control group object;
#'meaning samples where forward and reverse reads are counted separately 
#'cannot be in the same control group object as samples where forward and reverse
#'reads are counted together.
#'
#'A control group object with duplicate samples or samples with different 
#'correction statusses is possible but not recommended and will generate a 
#'warning message. 
#'
#'@return NIPTControlGroup object
#'
#'@examples 
#' \dontrun{
#' ##Retrieve filenames
#' bam_filepaths <- list.files(path = "/Path/to/bamfiles/", pattern = ".bam", full.names = T)
#' ##Load files and convert to control group
#' control_group  <- as_control_group(nipt_samples = lapply(X = bam_filepaths, bin_bam_sample, 
#'                                                          do_sort = F, separate_strands = FALSE))
#' ##Save control group for later
#' saveRDS(object = control_group, file = "/Path/to/directory/control_group.rds")
#' }
#'@export
as_control_group <- function(nipt_samples, control_group_type = generic_control_group){
#   if (length(unique(sapply(nipt_samples, getsamplenames))) != length(nipt_samples)){
#     cat("Warning, there appear to be duplicate sample names in control group \n")
#   }
   if ((length(unique(sapply(nipt_samples, getstrandtype))) != 1)){
      stop("More than one strand type in control group", call. = F)
   }
  control_group <- list()
  control_group[[samples]] <- nipt_samples
  control_group[[correction_status_autosomal_chromosomes]] <- unique(x = sapply(X = nipt_samples, FUN = getcorrectionstatus, 
                                                                            status_type = correction_status_autosomal_chromosomes))
  control_group[[correction_status_sex_chromosomes]] <- unique(x = sapply(nipt_samples, getcorrectionstatus,
                                                                      status_type = correction_status_sex_chromosomes))
  control_group[[description]] <- control_group_type
  
  class(control_group) <- c(NIPT_control_group_class, unique(sapply(nipt_samples, getstrandtype)))
  return(control_group)
}
#'Remove a sample by samplename from control group
#'
#'@param samplename Regular expression string. All matching samplenames 
#'are removed from the control group
#'@param nipt_control_group NIPTControlGroup object to remove samples from
#'@details
#'This function removes a sample from the `NIPTControlGroup` object by name.
#'Note that this function uses a regular expression, and if more sample_names 
#'satisfy the regular expression, they will also be removed. It returns a new 
#'NIPTControlGroup object. 
#'
#'@return NIPTControlGroup object
#'
#'@examples 
#' \dontrun{
#' new_control_group <- remove_sample_controlgroup(samplename = unwanted_sample, 
#'                                                 nipt_control_group = old_control_group)
#' }
#'
#'@export
remove_sample_controlgroup <- function(samplename, nipt_control_group){
  indices <- grep(pattern = samplename, x = sapply(nipt_control_group[[samples]], getsamplenames))
  if (length(indices != 0)){
  return (as_control_group(nipt_control_group[[samples]][-indices]))
  }
  else{
    message("No matching samplenames")
    return(nipt_control_group)
  }
}
#'Remove duplicate samples from control group
#'
#'Removes all duplicate samples in control group by samplename. 
#'@param nipt_control_group NIPTControlGroup object
#'
#'@details 
#'This functions removes duplicate samples from the control group based on name. 
#'It returns a new NIPTControlGroup object.
#'@return NIPTControlGroup object
#'
#'@examples 
#' \dontrun{
#' new_control_group <- remove_duplicates_controlgroup(nipt_control_group = old_control_group)
#' }
#'
#'@export
remove_duplicates_controlgroup <- function(nipt_control_group){
  indices <- which(duplicated(sapply(nipt_control_group[[samples]], getsamplenames)))
  as_control_group(nipt_control_group[[samples]][-indices])
}
#'Add a sample to an existing control group
#'
#'This functions adds NIPTSample objects to an existing control group and returns a new 
#'NIPTControlGroup object.
#'
#'@param nipt_control_group The NIPTControlGroup to add the samples to
#'@param samples_to_add A list with sample(s) to add. This always needs to be a list
#'
#'@return NIPTControlGroup object
#'
#'@examples 
#' \dontrun{
#' ##First bin the new sample
#' new_binned_sample <- bin_bam_sample(bam_filepath = "/path/to/file.bam", 
#'                                     separate_strands = T)
#' 
#' ##Then add the sample to the control group
#' new_control_group <- add_samples_controlgroup(nipt_control_group = my_control_group, 
#'                                               samples_to_add = new_binned_sample)
#' }
#'
#'@export
add_samples_controlgroup <- function(nipt_control_group, samples_to_add){
  as_control_group(nipt_samples = c(nipt_control_group[[samples]], samples_to_add))
}
