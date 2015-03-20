\name{probesetPermutations}
\alias{probesetPermutations}
\title{Permutations of the samples at each probeset}
\description{
    The function takes in the data.frame from 'firmaAnalysis' (containing 
    log2 FIRMA scores for all probe sets/exons (rows) in all samples 
    (columns)), along with the list 'percentiles' from 'aberrantExonUsage' 
    (containing the lower and upper percentile values of FIRMA scores used
    as thresholds for denoting aberrant exon usage), and makes 
    permutations of the FIRMA scores for each probe set/exon across all 
    samples. Based on the permutation, random relative amounts of aberrant 
    exon skipping and inclusion per sample is calculated and returned.}
\usage{
    probesetPermutations(fs, percentiles)
}
\arguments{
    \item{fs}{FIRMA scores for each probe set/sample combination (the 
        data.frame returned from the function 'firmaAnalysis').}
    \item{percentiles}{List containing two FIRMA score percentile values, 
        i.e., the lower and upper percentiles used as thresholds for 
        denoting aberrant exon usage (the list object 'percentiles' 
        returned from the function 'aberrantExonUsage').
    }
}
\value{
    A list with two vectors, with number of exon skipping and inclusion
    events, respectively, for each sample after permutations of the
    expression levels at each probeset.
}
\examples{
# Set up data set with FIRMA scores and calculate relative aberrant
# exon usage for each sample
    fs <- firmaAnalysis(useToyData=TRUE)
    tra <- aberrantExonUsage(1.0, fs)
# Make permutations of the expression data at each probeset
    perms <- probesetPermutations(fs, quantiles)
}
