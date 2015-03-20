\name{aberrantExonUsage}
\alias{aberrantExonUsage}
\title{Calculations of relative aberrant exon usage amounts per sample
    (based on FIRMA scores)}
\description{
    The function takes in the data.frame from 'firmaAnalysis' (containing 
    log2 FIRMA scores for all probe sets/exons (rows) in all samples 
    (columns)), and a number indicating which percentile value of global 
    FIRMA scores to be used as threshold for denoting aberrant exon 
    usage (default value '1', calculating the lower and upper 1st 
    percentiles, indicating aberrant exon skipping and inclusion, 
    respectively). Lower and upper percentile values are calculated and
    stored in the global list object 'quantiles'. Also, the total number 
    of exons per sample denoted with aberrant exon usage (having FIRMA 
    scores outside the indicated threshold values) is calculated and 
    stored in the global list object 'aberrantExons'. The function returns
    a vector with these total sample-wise amounts of aberrant exon usage
    (sum of aberrant skipping and inclusion amounts) relative to the
    average sample-wise amount in the dataset (log2-transformed).}
\usage{
    aberrantExonUsage(percentile, fs)
}
\arguments{
    \item{percentile}{This number indicates the percentile value of the 
        global FIRMA scores to be used as threshold for denoting aberrant 
        exon usage. Default value '1' calculates the lower and upper 1st
        percentiles, indicating aberrant exon skipping and inclusion, 
        respectively.}
    \item{fs}{Data.frame consisting of log2 FIRMA scores for all probe 
        sets/exons (rows) in all samples (columns). This data.frame is the 
        output from 'firmaAnalysis'.}
}
\value{
    A numeric vector with log2-transformed sample-wise amounts of aberrant 
    exon usage relative to the average sample-wise amount in the dataset.
    In addition, the quantiles list object is created, which contains the
    threshold values for the lower and upper percentiles.
}
\examples{
# Calculate aberrant exon usage for each sample in the data set:
    fs <- firmaAnalysis(useToyData=TRUE)
    tra <- aberrantExonUsage(1.0, fs)
#   The aberrantExonUsage function also creates the 'quantiles' object with
#   upper and lower threshold values for accepting aberrant exon usage, and
#   the list object 'aberrantExons' with the sample-wise number of exons
#   outside the threshold values.

}
