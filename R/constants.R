
library(Rsamtools)
library(sets)

Uncorrected <- "Uncorrected"
CombinedStrands <- "CombinedStrands"
SeparatedStrands <- 'SeparatedStrands'
GCcorrected <- "GC corrected"
ChiCorrected <- "Chi square corrected"
practical <- "Practical_CV"
theoretical <- "Theoretical_CV"
zscore <- "Z_score"
generic_control_group <- "General control group"
XY <- c("X", "Y")
cv <- "CV"
shapiro <- "Shapiro_P_value"
loess <- "LOESS"
bin <- "bin"

autosomal_chromosome_reads <- "autosomal_chromosome_reads"
sex_chromosome_reads <- "sex_chromosome_reads"
correction_status_autosomal_chromosomes <- "correction_status_autosomal_chromosomes"
correction_status_sex_chromosomes <- "correction_status_sex_chromosomes"
sample_name <- "sample_name"
samples <- "samples"
sample_names <- "sample_names"
description <- "description"


n_total_chromosomes <- 24
n_autosomal_chromosomes <- 22
autosomal_chromosomes <- 1:22
sex_chromosomes <- 23:24 
control_chromosomes <- as.character(c(1:12,14:17,19:20,22))
chromosomes_trisomy <- c(13,18,21)

n.models <- 4
n.predictors <- 4

bin_size <- 50000


rownames_combined_autosomal <- as.character(autosomal_chromosomes)
rownames_combined_sex <- XY
rownames_separated_forward_autosomal <- paste0(autosomal_chromosomes, "F")
rownames_separated_reverse_autosomal <- paste0(autosomal_chromosomes, "R")
rownames_separated_forward_sex <- paste0(XY, "F")
rownames_separated_reverse_sex <- paste0(XY, "R")

NIPT_sample_class <- "NIPTSample"
NCV_template_class <- "NCVTemplate"
NCV_result_class <- "NCVResult"
Z_template_class <- "ZscoreResult"
Regressions_result_class <- "RegressionResult"
NIPT_control_group_class <- "NIPTControlGroup"



