\name{scatterPlot}
\alias{scatterPlot}
\title{Scatterplot showing relative amounts of aberrant exon usage per 
    sample
}
\description{
    Scatterplot showing sample-wise relative amounts (blue dots) of 
    aberrant exon inclusion (horizontal axis) and exon skipping (vertical 
    axis). Random sample-wise amounts calculated from permuted FIRMA 
    scores can also be included in the plot (yellow dots).}
\usage{
    scatterPlot(fileName, permutations, percentileHits, 
    permPercentileHits)
}
\arguments{
    \item{fileName}{Output filename. File format is optional, but must be
        one of png, jpg, eps or pdf.}  
    \item{permutations}{Boolean argument to indicate whether permuted data
        should be included. 'TRUE' adds permutation data to the plot.}
    \item{percentileHits}{List with two items containing sample-wise 
        numbers of exons denoted with aberrant exon skipping and inclusion, 
        i.e., having FIRMA scores in the indicated lower and upper 
        percentiles, respectively. This list is returned from the function
        'aberrantExonUsage'.}
    \item{permPercentileHits}{List with two items containing random 
        sample-wise numbers of exons denoted with aberrant exon skipping and 
        inclusion, i.e., having FIRMA scores in the indicated lower and upper 
        percentiles, respecively (calculated from FIRMA score permutations).
        This list is returned from the function 'probesetPermutations'.}
    }
\value{
    scatterPlot is used for the side-effect of producing a scatter plot
    showing relative amounts of aberrant exon usage per sample.
}
\examples{
    data(splicingFactors)
    fs <- firmaAnalysis(useToyData=TRUE)
    gs <- readGeneSummaries(useToyData=TRUE)
    tra <- aberrantExonUsage(1.0, fs)
#   The aberrantExonUsage function also creates the 'quantiles' object with
#   upper and lower threshold values for accepting aberrant exon usage, and
#   the list object 'aberrantExons' with the sample-wise number of exons
#   outside the threshold values.
    aberrantExonsPerms <- probesetPermutations(fs, quantiles)

# Create scatter plot with the samples
    scatterPlot("scatter.png", TRUE, aberrantExons, aberrantExonsPerms)

}
