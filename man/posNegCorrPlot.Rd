\name{posNegCorrPlot}
\alias{posNegCorrPlot}
\title{posNegCorrPlot}
\description{
    The posNegCorrPlot is a scatterPlot that compares the amount of splicing 
    factor genes (red) for which expression levels are significant 
    positively (vertical axis) and negatively (horizontal axis) correlated 
    with the total relative amounts of aberrant exon usage per sample. The 
    plot can also include results from permutations of the sample-wise 
    aberrant exon usage amounts (dark blue), and randomly constructed gene 
    sets of 280 genes.
}
\usage{
    posNegCorrPlot(fileName, tra, geneSummaries, splicingFactors, 
        randomGeneSets, traPermutations)
}
\arguments{
    \item{fileName}{Output filename. File format is optional, but must be
        one of png, jpg, eps or pdf.}
    \item{tra}{List object with total relative amounts of aberrant exon 
        usage per sample. The object is returned by the function 
        'aberrantExonUsage'.}
    \item{geneSummaries}{The data.frame with gene-level expression values 
        for each sample, returned from the function 'readGeneSummaries'.}
    \item{splicingFactors}{List of genes (Affymetrix transcript 
        cluster id's and gene symbols) involved in pre-mRNA splicing. An
        example set of 280 genes is obtained by issuing the command 
        'data(splicingFactors)', but the input list can include any set of
        genes included in the data.frame returned from the function
        'readGeneSummaries'.}
    \item{randomGeneSets}{Number of random gene sets of 280 genes to be 
        created and included in the analysis.}
    \item{traPermutations}{Number of permutations of the sample-wise
        amounts of aberrant exon usage to be performed and included in the 
        analysis.}
}
\value{
    posNegCorrPlot is used for the side-effect of producing a scatter plot
    that compares the amount of splicing factor genes (red) for which
    expression levels are significant positively (vertical axis) and
    negatively (horizontal axis) correlated with the total relative amounts
    of aberrant exon usage per sample. In addition, the plot can also
    include results from permutations of the sample-wise aberrant exon usage
    amounts (dark blue), and randomly constructed gene sets of 280 genes.
}
\examples{
    data(splicingFactors)
    fs <- firmaAnalysis(useToyData=TRUE)
    gs <- readGeneSummaries(useToyData=TRUE)
    tra <- aberrantExonUsage(1.0, fs)

# Create plot that compares the amount of splicing factor genes for which
# expression levels are significant positively and negatively correlated
# with the total relative amounts of aberrant exon usage per sample
    posNegCorrPlot("cg.png", tra, gs, splicingFactors, 1000, 1000)
}
