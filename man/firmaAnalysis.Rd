\name{firmaAnalysis}
\alias{firmaAnalysis}
\title{Read CEL files and perform FIRMA analysis}
\description{
    The function makes use of the aroma.affymetrix package to analyze 
    Affymetrix Human Exon 1.0 ST Arrays. The function reads CEL files, and 
    performs background correction, normalization (customized RMA 
    approach), and alternative splicing analysis according to the FIRMA 
    method 
    (http://www.aroma-project.org/vignettes/FIRMA-HumanExonArrayAnalysis).  
    The function returns a data.frame with log2 FIRMA (alternative 
    splicing) scores for each probeset/sample combination.
}
\usage{
    firmaAnalysis(useToyData, aromaPath, dataSetName)
}
\arguments{
    \item{useToyData}{Boolean argument to indicate whether sample data sets
        included in the TIN package should be used in the analysis.}
    \item{aromaPath}{Absolute or relative path to the aroma.affymetrix
        directory. Requires custom CDF annotation file (please refer to the
        FIRMA vignette for download and setup).}
    \item{dataSetName}{Name of folder in the 'aromaPath' containing raw
        data (CEL files; please refer to the FIRMA vignette for setup).}
}
\value{
    A data.frame with expression level values after the FIRMA analysis has been
    applied. The data.frame consists of one column for each sample and one row 
    for each probeset.
}
\references{
    E. Purdom, K. Simpson, M. Robinson, J. Conboy, A. Lapuk & T.P. Speed, 
    FIRMA: a method for detection of alternative splicing from exon array 
    data. Bioinformatics, 2008.
}
\examples{
# Perform FIRMA analysis on the raw expression data
# To use sample data sets included in the TIN package as input:
    fs <- firmaAnalysis(useToyData=TRUE)

# To use your own data, provide path to aroma.affymetrix root directory and
# name of data set as arguments:
#   fs <- firmaAnalysis(useToyData=FALSE, "/tmp/path/to/aroma.affymetrix",
#   "sampleSet")

}

