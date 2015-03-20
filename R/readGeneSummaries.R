####################################################################
## Author: Bjarne Johannessen, Anita Sveen and Rolf I. Skotheim
## Maintainer: Bjarne Johannessen <bjarnej@ifi.uio.no>
## License: Artistic 2.0
## Part of the TIN package
####################################################################

## Function for reading pre-pocessed and summarized gene-level expression
## data from file, and building a data.frame with the information.

readGeneSummaries <- function(useToyData=FALSE, summaryFile)
{
    if (useToyData) {
        sampleSetGeneSummaries <- NULL
        data(sampleSetGeneSummaries, envir=environment())
        geneSummaries <- sampleSetGeneSummaries
    } else {
        if (missing(summaryFile)) {
            stop("Please provide file with gene level expression data",
                call.=FALSE)
        } else {
        geneSummaries <- 
            read.table(summaryFile, sep="\t", header=TRUE, row.names=1)
        }
    }
    return(geneSummaries)
}
