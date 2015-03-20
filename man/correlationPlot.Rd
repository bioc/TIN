\name{correlationPlot}
\alias{correlationPlot}
\title{correlationPlot}
\description{
    Function for creating a plot that visualizes the number of splicing
    factor genes with expression levels significantly correlated with the
    sample-wise total relative amounts of aberrant exon usage (red). To
    compare this value with correlation values obtained from random sample
    permutations, the function performs two types of iterative sample
    calculations as control experiments. First, expression levels of
    the splicing factor gene set are correlated with permutations of the
    relative aberrant exon usage (dark blue). Second, expression levels of
    randomly generated gene sets are correlated with the relative aberrant
    exon usage amounts in the data (pale blue).
}
\usage{
    correlationPlot(fileName, tra, geneSummaries, splicingFactors, 
    randomGeneSets, traPermutations)
}
\arguments{
    \item{fileName}{Output filename. File format is optional, but must be 
        one of png, jpg, eps or pdf.}
    \item{tra}{List of sample-wise total relative amounts of aberrant exon
        usage, obtained using the 'aberrantExonUsage' function.}
    \item{geneSummaries}{Data.frame with gene-level expression data for 
        each sample, returned from the function 'readGeneSummaries'.}
    \item{splicingFactors}{List with Affymetrix transcript cluster id's 
        and gene symbols for a set of genes involved in pre-mRNA splicing. An 
        example set with 280 genes is obtained by issuing the command
        'data(splicingFactors)', but the list may well include any set of 
        genes included in the data.frame returned from the function
        'readGeneSummaries'.}
    \item{randomGeneSets}{Number of random gene sets of 280 genes to be 
        created and included in the analysis.}
    \item{traPermutations}{Number of permutations of the sample-wise 
        amounts of aberrant exon usage to be performed and included in the 
        analysis.}
}
\value{
    correlationPlot is used for the side-effect of producing a plot showing
    the number of splicing factor genes with expression levels significantly 
    correlated with the sample-wise total relative amounts of aberrant exon
    usage.
}
\examples{
    data(splicingFactors)
    fs <- firmaAnalysis(useToyData=TRUE)
    gs <- readGeneSummaries(useToyData=TRUE)
    tra <- aberrantExonUsage(1.0, fs)

# Create a plot that visualizes the number of splicing factor
# genes with expression levels significantly correlated with the
# sample-wise total relative amounts of aberrant exon usage.
    correlationPlot("c.png", tra, gs, splicingFactors, 1000, 1000)

}
