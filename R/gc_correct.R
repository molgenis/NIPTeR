#'Perform a GC bias correction on nipt sample
#'
#'LOESS based GC bias correction algorithm described by Chen et al (2011)
#'
#'@param nipt_object The object that will be corrected. This can either be a 
#'`NIPTSample` or a `NIPTControlGroup` object
#'@param method To select the LOESS based method use  \emph{"LOESS"},  
#'to select the bin weights based method use  \emph{"bin"}. 
#'@param include_XY Also apply correction to X and Y chromosomes? 
#'@param span The span for the LOESS fit. Only applicable when LOESS method is used. 
#'@param ref_genome The reference genome used. Either \emph{"hg37"} or \emph{"hg38"}
#'default = \emph{"hg37"}
#'@details
#'GC content bias is the correlation between the number of reads mapped to a specific genomic 
#'region and the GC content of this region. In NIPTeR, two GC bias correction algorithms 
#'have been implemented, the LOESS based method introduced by Chen et al. (2011) and the bin 
#'weight based method described by Fan and Quake (2010). 
#'@return Depending on the input object either a NIPTSample or a NIPTControlGroup object 
#'
#'@examples 
#' \dontrun{
#' ##Correct NIPTSample object using LOESS method
#' loess_corrected_sample <- gc_correct(nipt_object = sample_of_interest, method = "LOESS",
#'                                      include_XY = F, span = 0.75)
#' ##Correct NIPTControlGroup object using bin method
#' gc_bin_corrected_control_group <- gc_correct(nipt_object = control_group, method = "bin", 
#'                                              include_XY = T)
#' }
#' 
#'@export
gc_correct <- function(nipt_object, method = "LOESS", include_XY = F, span =0.75, 
                       ref_genome = "hg37"){
  if(class(nipt_object)[1] == NIPT_sample_class){
    if (method == loess){
       return(corrected_sample <- gc_correct_NIPTSample_loess(nipt_object, span = 0.75, 
                                                              include_XY = include_XY, 
                                                              ref_genome = ref_genome))
    }
    if (method == bin){
      
      return(gc_correct_NIPTSample_bin(nipt_object, 
                                       include_XY = include_XY, ref_genome = ref_genome))
    }
    else{
      stop("Error")
    }
  }
  if(class(nipt_object)[1] == NIPT_control_group_class){
    if (method == loess){
      return(gc_correct_NIPTControlGroup_loess(nipt_object, span = 0.75, include_XY = include_XY,
                                               ref_genome = ref_genome))
    }
    if (method == bin){
      return(gc_correct_NIPTControlGroup_bin(nipt_object, include_XY = include_XY,
                                             ref_genome = ref_genome))
    }
    else{
      stop("Error")
    }
  }
}