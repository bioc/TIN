\name{readGeneSummaries}
\alias{readGeneSummaries}
\title{Read gene-level expression summaries}
\description{
    The function reads pre-pocessed and summarized gene-level expression 
    data from file, and builds a data.frame with the information. The 
    data.frame will have the same structure as the input summary file, 
    i.e., genes in rows and samples in columns.}
\usage{
    readGeneSummaries(useToyData, summaryFile)
}
\arguments{
    \item{useToyData}{Boolean argument to indicate whether sample data sets
        included in the TIN package should be used in the analysis.}
    \item{summaryFile}{Tab separated file with expression values for each 
        combination of gene (row) and sample (column). Gene expression 
        summary files can be obtained by using for instance APT (Affymetrix
        Power Tools) or Expression Console. The file is expected to have one
        initial header line with sample names, and one initial column with
        Affymetrix transcript cluster id's as row names.
    }
}
\value{
    A data.frame with gene summary values. It consists of one column for each
    sample, and one row for each gene.
}
\examples{
# Read pre-processed gene summary values from file
# To use sample test data included in the TIN package as input:
    gs <- readGeneSummaries(useToyData=TRUE)
# To use your own data, provide path to txt file with expression values:
#   gs <- readGeneSummaries("/tmp/path/to/GeneLevelExpressionValues.txt")
}
