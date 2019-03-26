## Update v1.0.2 2016-03-09
This is a minor update: In this version I have:

Changed, due to a change in S4Vectors, removed "#'@importFrom stats sd" from bin_sample.R, so that stats::sd will not replace the 
previous import done by S4Vectors (which created a warning). 
In the following packages "sd" has been replaced by "stats::sd":
diagnose_control_group.R
prepare_ncv_template.R
regression.R
z_score.R

## Update v1.0.1 2015-12-29
This is a minor update: In this version I have:

Added contact information to vignette
Added examples to the function information fields (\dontrun{}, since large datasets would be necessary to test)

## Resubmission
This is a resubmission (2). In this version I have:

Changed title from "Fast and Accurate Trisomy Prediction in NIPT" 
to "Fast and Accurate Trisomy Prediction in Non-Invasive Prenatal
    Testing".

This is a resubmission. In this version I have:

* Changed title from "R Package for Non-Invasive Prenatal Testing Analysis"" to NIPTeR: Fast and Accurate Trisomy Prediction in NIPT.

* Fixed errors "Undefined global functions or variables:
  lm loess median predict sd shapiro.test" by importing the functions
  to the namespace

##  Test environments
* OS X x86_64-apple-darwin13.4.0 (64-bit) | R 3.2.2
* Windows 7 Enterprise SP 1 (32-bit) | R 3.2.2
* devtools::build_win()


## R CMD check results
There were no ERRORs, WARNINGs or NOTES. 

##Update Sosa
Added to constants.R 
bin_sequence <-
  structure(c(20L, 27L, 28L, 30L, 33L, 37L, 42L, 44L, 47L, 1L, 
              2L, 4L, 5L, 6L, 7L, 8L, 10L, 15L, 17L, 23L, 24L, 26L, 92L, 93L
  ), .Names = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", 
                "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", 
                "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", 
                "chr22", "chrX", "chrY"))

Changed splitted reads in bin_sample.R
  splitted_reads[[1]] <- splitted_reads[[1]][bin_sequence]
  splitted_reads[[2]] <- splitted_reads[[2]][bin_sequence]