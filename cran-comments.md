##Update v1.0.3 2019-04-01
This is a minor update: In this version I have:
changed the bin_sample.R script to be able to process lexicographically sorted bam files.


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

