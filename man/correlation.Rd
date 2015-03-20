\name{correlation}
\alias{correlation}
\title{Calculates the correlation between sample-wise amounts of aberrant 
    exon usage and splicing factor expression levels} 
\description{
    The function makes use of the corAndPValue function from the WGCNA
    package to calculate sample-wise Pearson correlation between relative
    amounts of aberrant exon usage and splicing factor expression 
    levels.}
\usage{
    correlation(splicingFactors, geneSummaries, tra)
}
\arguments{
    \item{splicingFactors}{A data.frame with a list of splicing factor  
        genes (Affymetrix transcript cluster id's and gene symbols)
        to be included in the correlation analysis. The list can 
        include any set of genes included in the data.frame returned from 
        the function 'readGeneSummaries'. An example set of 280 genes is
        obtained by issuing the command 'data(splicingFactors)'.}
    \item{geneSummaries}{A data.frame with gene-level expression values 
        for all samples, returned from the 'readGeneSummaries' function.}
    \item{tra}{List of sample-wise total relative amounts of aberrant exon 
        usage, returned from the 'aberrantExonUsage' function.
    }
}
\value{
    A list with sample-wise Pearson correlation values between relative
    amounts of aberrant exon usage and splicing factor expression
    levels.
}
\examples{
    data(sampleSetGeneSummaries)
    gs <- sampleSetGeneSummaries
    fs <- firmaAnalysis(useToyData=TRUE)
    tra <- aberrantExonUsage(1.0, fs)
# calculate correlation between splicing factor expression and aberrant
# exon usage
    data(splicingFactors)
    corr <- correlation(splicingFactors, gs, tra)
}
