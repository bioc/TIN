test_aberrantExonUsage <- function() {
    checkEquals(dim(firmaAnalysis(useToyData=TRUE))[2], length(tra <- aberrantExonUsage(1.0, firmaAnalysis(useToyData=TRUE))))
    checkTrue(!anyNA(aberrantExonUsage(1.0, firmaAnalysis(useToyData=TRUE))))
    checkTrue(as.numeric(quantiles[1]) < 0)
    checkTrue(as.numeric(quantiles[2]) > 0)
}

