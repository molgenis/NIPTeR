% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/control_group.R
\name{add_samples_controlgroup}
\alias{add_samples_controlgroup}
\title{Add a sample to an existing control group}
\usage{
add_samples_controlgroup(nipt_control_group, samples_to_add)
}
\arguments{
\item{nipt_control_group}{The NIPTControlGroup to add the samples to}

\item{samples_to_add}{A list with sample(s) to add. This always needs to be a list}
}
\value{
NIPTControlGroup object
}
\description{
This functions adds NIPTSample objects to an existing control group and returns a new 
NIPTControlGroup object.
}
\examples{
\dontrun{
##First bin the new sample
new_binned_sample <- bin_bam_sample(bam_filepath = "/path/to/file.bam", 
                                    separate_strands = T)

##Then add the sample to the control group
new_control_group <- add_samples_controlgroup(nipt_control_group = my_control_group, 
                                              samples_to_add = new_binned_sample)
}

}

