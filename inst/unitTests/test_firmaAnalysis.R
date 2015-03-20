test_firmaAnalysis <- function() {
    checkEquals(sum(sapply(firmaAnalysis(useToyData=TRUE), is.numeric) == TRUE), ncol(firmaAnalysis(useToyData=TRUE)))
    checkTrue(!anyNA(firmaAnalysis(useToyData=TRUE)))
}

