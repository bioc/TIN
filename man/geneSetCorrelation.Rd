\name{geneSetCorrelation}
\alias{geneSetCorrelation}
\title{Correlation between aberrant exon usage and expression levels for a
    number of gene sets.
}
\description{
    The function makes use of the corAndPValue function from the WGCNA 
    package to calculate the Pearson correlation between sample-wise
    aberrant exon usage amounts and expression levels of all genes for all
    gene sets defined by the input parameter list geneSets.}
\usage{
    geneSetCorrelation(geneSets, geneAnnotation, geneSummaries, tra, 
        noGeneSets)
}
\arguments{
    \item{geneSets}{A data.frame with a number of gene lists to be included
        in the correlation analysis. An example data.frame of 1,454 lists of
        genes (specified by Affymetrix transcript cluster id's and gene
        symbols) is included in the package, and will be accessible by
        issuing the command 'data(genesets)'. The example data.frame contains
        a complete collection of all Gene Ontology gene sets included in the
        Molecular Signatures Database v3.1.}
    \item{geneAnnotation}{A data.frame with Affymetrix transcript cluster 
        id's and gene symbols for all 22,011 genes included in the 
        Affymetrix 'core' set.}
    \item{geneSummaries}{A data.frame with gene-level expression values
        for all samples, returned from the 'readGeneSummaries' function.}
    \item{tra}{List with total relative amounts of aberrant exon usage for all
        samples obtained using the 'aberrantExonUsage' function.}
    \item{noGeneSets}{Optional argument specifying how many gene sets to
        include in the analysis.}
}
\value{
    A data.frame with one row for each data set used as input, and columns
    for name of set, number of genes, number of significant
    positively/negatively correlated genes in the set, and median correlation
    strength.
}
\examples{
# Load data
    data(geneSets)
    data(geneAnnotation)
    fs <- firmaAnalysis(useToyData=TRUE)
    gs <- readGeneSummaries(useToyData=TRUE)
    tra <- aberrantExonUsage(1.0, fs)
# Calculate correlation in other gene sets
    crs <- geneSetCorrelation(geneSets, geneAnnotation, gs, tra, 50)
}
