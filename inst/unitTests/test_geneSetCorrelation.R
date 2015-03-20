test_geneSetCorrelation <- function() {
    data(geneSets)
    data(geneAnnotation)
    checkTrue(!anyNA(geneSetCorrelation(geneSets, geneAnnotation, readGeneSummaries(useToyData=TRUE), aberrantExonUsage(1.0, firmaAnalysis(useToyData=TRUE)), 20)[, 3]))
    checkTrue(!anyNA(geneSetCorrelation(geneSets, geneAnnotation, readGeneSummaries(useToyData=TRUE), aberrantExonUsage(1.0, firmaAnalysis(useToyData=TRUE)), 20)[, 4]))
    checkTrue(!anyNA(geneSetCorrelation(geneSets, geneAnnotation, readGeneSummaries(useToyData=TRUE), aberrantExonUsage(1.0, firmaAnalysis(useToyData=TRUE)), 20)[, 5]))
    checkEquals(sum(sapply(geneSetCorrelation(geneSets, geneAnnotation, readGeneSummaries(useToyData=TRUE), aberrantExonUsage(1.0, firmaAnalysis(useToyData=TRUE)), 20), is.numeric) == TRUE), 4)
    checkTrue(all(geneSetCorrelation(geneSets, geneAnnotation, readGeneSummaries(useToyData=TRUE), aberrantExonUsage(1.0, firmaAnalysis(useToyData=TRUE)), 20)[, 5] > -1))
    checkTrue(all(geneSetCorrelation(geneSets, geneAnnotation, readGeneSummaries(useToyData=TRUE), aberrantExonUsage(1.0, firmaAnalysis(useToyData=TRUE)), 20)[, 5] < 1))
}

