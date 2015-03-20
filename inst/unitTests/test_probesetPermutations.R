test_probesetPermutations <- function() {
    aberrantExonUsage(1.0, firmaAnalysis(useToyData=TRUE))
    perms <- probesetPermutations(firmaAnalysis(useToyData=TRUE), quantiles)
    checkEquals(sum(as.numeric(perms[[1]]))/100, ncol(readGeneSummaries(useToyData=TRUE)))
    checkEquals(sum(as.numeric(perms[[1]]))/100, sum(sapply(firmaAnalysis(useToyData=TRUE), is.numeric) == TRUE))
    checkEquals(sum(as.numeric(perms[[2]]))/100, ncol(readGeneSummaries(useToyData=TRUE)))
    checkEquals(sum(as.numeric(perms[[2]]))/100, sum(sapply(firmaAnalysis(useToyData=TRUE), is.numeric) == TRUE))
    checkTrue(sum(as.numeric(perms[[1]])) %% 100 == 0)
    checkTrue(sum(as.numeric(perms[[2]])) %% 100 == 0)
}

