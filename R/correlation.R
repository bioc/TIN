####################################################################
## Author: Bjarne Johannessen, Anita Sveen and Rolf I. Skotheim
## Maintainer: Bjarne Johannessen <bjarnej@ifi.uio.no>
## License: Artistic 2.0
## Part of the TIN package
####################################################################

## Function to calculate the correlation between sample-wise amounts of
## aberrant exon usage and splicing factor expression levels.



correlation <- function(splicingFactors, geneSummaries, tra)
{

## Select expression values from 280 splicing factors.
    spliceGeneSummaries <- 
        geneSummaries[rownames(geneSummaries) %in% splicingFactors[, 1], ]

    correlation <- corAndPvalue(t(spliceGeneSummaries[, seq_along(tra)]), 
        tra, alternative = "two.sided")

    return(correlation)

}
