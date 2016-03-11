.unlist <- function (x){
  x1 <- x[[1L]]
  if (is.factor(x1)){
    structure(unlist(x), class = "factor", levels = levels(x1))
  } else {
    do.call(c, x)
  }
}
#' Load and bin BAM file 
#' 
#' Load a BAM file and count reads in bins of size 50.000 base pairs
#'@import sets
#'@import Rsamtools
#'@import S4Vectors
#'@importFrom stats lm
#'@importFrom stats loess
#'@importFrom stats median
#'@importFrom stats predict
#'@importFrom stats shapiro.test
#'@param bam_filepath Character The location and filename on the file system where the bam file is stored
#'@param do_sort Boolean Sort the bam file? If the bam is unsorted set to true, 
#'but the use of pre-sorted bam files is recommended.
#'@param separate_strands Boolean If set to true, reads from forward and reverse strands are counted and stored separately. 
#'This option should be used if you are planning on using regression, since this doubles the number of 
#'predictors (F+R) and distributes predictive power more equally over prediction sets since F and R 
#'strand from the same chromosome cannot be both in one prediction set. 
#'@param custom_name String The name of sample. Default samplename is the filename of the bam 
#'file without the .bam suffix and filepath prefix. 
#'@aliases nipt_sample
#'
#'@return Object NIPTSample
#'
#'@details
#'This function returns an object of class NIPTSample, the main 'currency' of this package.
#' It is a list with 5 items:
#'\itemize{
#' \item List \strong{autosomal_chromosome_reads} Autosomal reads are stored in
#' a matrix where the columns are the bins and rows (22) represent the autosomal
#' chromosomes. The length of this list is either 1 or 2, depending if the 
#' forward and reverse reads are counted separately.
#' \item Character \strong{correction_status_autosomal_chromosomes} The correction
#' status of the autosomal reads. The status can either be \emph{Uncorrected} or
#' \emph{GC Corrected} and/or \emph{Chi Corrected}
#' \item List \strong{sex_chromosome_reads} Sex chromosome reads are stored in a 
#' similar matrix(es) as the autosomal chromosome reads, now with 2 (X and Y) rows.
#' \item Character \strong{correction_status_autosomal_chromosomes} The status can 
#' either be \emph{Uncorrected} or \emph{GC Corrected} and/or \emph{Chi Corrected}.
#' \item Character \strong{sample_name} Sample name 
#'}
#'
#'@examples 
#' \dontrun{
#' ##To process a single sample
#' binned_sample <- bin_bam_sample(bam_filepath = "/path/to/file.bam", 
#'                                 separate_strands = T)
#' 
#' ##To create a control group out of a set of bam files
#' bam_filepaths <- list.files(path = "/Path/to/bamfiles/", 
#'                             pattern = ".bam", full.names = T)
#' 
#' control_group  <- as_control_group(nipt_samples = lapply(X = bam_filepaths, 
#'                                    bin_bam_sample, do_sort = F, 
#'                                    separate_strands = T))
#' }
#' 
#'@export
bin_bam_sample <- function(bam_filepath, do_sort=FALSE, separate_strands=FALSE, custom_name = NULL){
  message("Loading Bam")
  if (!is.null(custom_name)){
    name = custom_name
  }
  else{
    name = basename(sub("^([^.]*).*", "\\1", bam_filepath))
  }
  if (do_sort == TRUE){
    temp <- sortBam(file = bam_filepath, destination = tempfile())
    bam <- scanBam(file = temp)
  }
  else{
    bam <- scanBam(bam_filepath)
  }
  message("BAM loaded")
  
  chromos <- bam[[1]][3][[1]]
  strands <- bam[[1]][4][[1]]
  reads <- bam[[1]][5][[1]]
  splitted_reads <- split(x = reads, f =  droplevels(strands[strands != "*"]))
  splitted_chromos <- split(x = chromos, f = droplevels(strands[strands != "*"]))
  binned_reads <- list()
  message("Binning")
  binned_reads[[1]] <- bin_reads(reads_data_frame = splitted_reads[[1]], chroms = splitted_chromos[[1]])
  binned_reads[[2]] <- bin_reads(reads_data_frame = splitted_reads[[2]], chroms = splitted_chromos[[2]])
  message("Binning done")
  
  if (separate_strands == FALSE){
    binned_reads <- list(Reduce("+", binned_reads))
  }
  
  autosomal_reads <- lapply(X = binned_reads, FUN = splitchromosomes, chromosomes = autosomal_chromosomes)
  sex_reads <- lapply(X = binned_reads, FUN = splitchromosomes, chromosomes = sex_chromosomes)
  new_sample <- construct_sample(autosomal_reads = autosomal_reads, sex_reads = sex_reads, name = name, 
                                 correction_status_autosomal = Uncorrected, correction_status_sex = Uncorrected)
  
  return(new_sample)
}

bin_reads <- function(reads_data_frame, chroms){
  n_bins <- getbins(max(reads_data_frame), bin_size)
  bin <- matrix(data = 0, nrow = 0, ncol = n_bins , dimnames = list(NULL, 1:n_bins))
  if ((length(unique(chroms[1:100]))) != 1){
    stop("BAM file appears to be unsorted", call. = F)
  }
  chromos <- rle(x = as.numeric(chroms))
  min_read <- 0
  max_read <- 1
  for (chromo in 1:n_total_chromosomes){
    max_read <- sum(chromos$lengths[1:(chromo)]) 
    
    reads <- sapply(X = unique(reads_data_frame[min_read:max_read]),  
                    FUN =  getbins, bin_size=bin_size)
    min_read <- max_read +1
    bins <- tabulate(reads, nbins = n_bins)
    bin <- rbind(bin, bins)
  }
  return (bin)
}
getbins <- function(pos, bin_size){
  (pos-1) %/% bin_size + 1
}