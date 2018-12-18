INTRODUCTION
------------

Non-Invasive Prenatal Testing (NIPT) is a method based on analysis of cell-free fetal (cff) DNA in maternal blood using next-generation sequencing. The basic approach is relatively straightforward: cell-free DNA is isolated from the mothers blood plasma and the obtained DNA fragments are sequenced. Subsequently the number of DNA fragments originating from the different chromosomes is counted, since in case of a fetal trisomy a relative increase of the fraction of reads of the trisomic chromosome is expected.

Various analysis methods for trisomy prediction, variation reduction and quality control have been described in the literature. These methods can be combined to improve the sensitivity of the analysis. The NIPTeR package makes these methods available and allows for flexibility in workflow design.

NIPTeR uses BAM files as input and calculates chromosomal fractions, which are used for trisomy prediction. These fractions are based on read counts of the different chromosomes. The chromosomes are divided in bins to allow for read count correction.

INSTALLATION INSTRUCTIONS
-------------------------

NIPTeR is dependent on three packages: sets, RSamtools and S4Vectors. For installation use:

install.packages("sets") source("<http://bioconductor.org/biocLite.R>") biocLite("Rsamtools") install.packages("NIPTeR")

EXAMPLES AND OVERVIEW
---------------------

Workflow examples and components overview are included in the NIPTeR vignette.
A step-by-step use case is added to the NIPTeR package paper as Additional file 1. This case makes use of the test samples provided in this git repository.

Literature
---------------------
NIPTeR methods are described in our [paper](https://www.nature.com/articles/s41598-017-02031-5).

Johansson L.F, de Boer E.N., de Weerd H.A., van Dijk F., Elferink M.G., Schuring-Blom G.H., Suijkerbuijk R.F., Sinke R.J., te Meerman G.J., Sijmons R.H., Swertz M.A. and Sikkema-Raddatz B. “Novel Algorithms for Improved Sensitivity in Non-Invasive Prenatal Testing,” Sci. Rep., vol. 7, no. 1, p. 1838, 2017.

Implementation of the NIPTeR package is described in our [paper](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2557-8)

Johansson L.F., de Weerd H.A., de Boer E.N., van Dijk F., te Meerman G.J., Sijmons R.H., Sikkema-Raddatz B. and Swertz M.A. , “NIPTeR : an R package for fast and accurate trisomy prediction in non-invasive prenatal testing,” BMC Bioinformatics, vol. 19, no. 531, pp. 1–5, 2018.
