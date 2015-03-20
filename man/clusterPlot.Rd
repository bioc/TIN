\name{clusterPlot}
\alias{clusterPlot}
\title{clusterPlot}
\description{
    Create plot from hierarchical clustering analysis of the samples, 
    based on splicing factor expression levels.}
\usage{
    clusterPlot(geneSummaries, tra, distmethod, clustermethod, fileName)
}
\arguments{
    \item{geneSummaries}{The data.frame with gene-level expression values 
        for each sample, returned from the function 'readGeneSummaries'.}
    \item{tra}{The list returned from the function 'aberrantExonUsage', 
        containing sample-wise total relative amounts of aberrant exon 
        usage.}
    \item{distmethod}{Which distance measure to be used. Possible options 
        are "euclidean", "maximum", "manhattan", "canberra", "binary" or 
        "minkowski".}
    \item{clustermethod}{Which clustering algorithm to be used. Possible 
        options are "ward", "single", "complete", "average", "mcquitty", 
        "median" or "centroid".}
    \item{fileName}{Output filename. File format is optional, but must be 
        one of png, jpg, eps or pdf.}
}
\value{
    clusterPlot is used for the side-effect of producing a hierarchical
    clustering plot showing how the samples are separated based on
    expression levels for the splicing factors in each sample.
}
\examples{
    fs <- firmaAnalysis(useToyData=TRUE)
    gs <- readGeneSummaries(useToyData=TRUE)
    tra <- aberrantExonUsage(1.0, fs)
# create cluster plot with the samples
    clusterPlot(gs, tra, "euclidean", "complete", "cluster.png")
}
