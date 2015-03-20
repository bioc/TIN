test_readGeneSummaries <- function() {
    checkEquals(sum(sapply(firmaAnalysis(useToyData=TRUE), is.numeric) == TRUE), ncol(readGeneSummaries(useToyData=TRUE)))
    checkTrue(!anyNA(readGeneSummaries(useToyData=TRUE)))
}

