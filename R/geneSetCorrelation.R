####################################################################
## Author: Bjarne Johannessen, Anita Sveen and Rolf I. Skotheim
## Maintainer: Bjarne Johannessen <bjarnej@ifi.uio.no>
## License: Artistic 2.0
## Part of the TIN package
####################################################################

## Function for calculating correlation between aberrant exon usage 
## and expression levels for a number of gene sets.




geneSetCorrelation <- function(geneSets, geneAnnotation, geneSummaries,
    tra, noGeneSets)
{
    if(missing(noGeneSets) || noGeneSets > dim(geneSets)[1]) {
        noGeneSets <- dim(geneSets)[1]
    }
    geneSetNames <- vector('character', noGeneSets)
    sizeGeneSets = vector('numeric', noGeneSets)
    noOfSignCorrGenes = vector('numeric', noGeneSets)
    noOfSignNegCorrGenes = vector('numeric', noGeneSets)
    medianCorrelationStrength <- vector('numeric', noGeneSets)
    for (i in seq_len(noGeneSets)) {
        geneSet <- unique(geneSets[i, ][geneSets[i, ] != ""])
        geneSetNames[i] <- geneSet[1]
        geneSet <- geneSet[3:length(geneSet)]
        transcriptClusterIds <- 
            geneAnnotation[geneAnnotation$Gene_symbol %in% geneSet, 1]

        sizeGeneSets[i] <- as.numeric(length(transcriptClusterIds))
        selGeneSummaries <- 
            geneSummaries[rownames(geneSummaries) %in% transcriptClusterIds, ]
        correlation <- corAndPvalue(t(selGeneSummaries[, seq_along(tra)]), 
            tra, alternative = "two.sided")
        m<-matrix(c(correlation$cor, correlation$p), ncol = 2)
        colnames(m) <- c("cor", "p")
        m.df <- data.frame(m)
        medianCorrelationStrength[i] <- median(m.df$cor)
        s <- m.df[m.df$p < 0.05, ]
        noOfSignCorrGenes[i] <- length(s$cor[s$cor > 0])
        noOfSignNegCorrGenes[i] <- length(s$cor[s$cor < 0])
    }
    A <- as.data.frame(matrix(nrow = noGeneSets, ncol = 5))
    colnames(A) <- c("Geneset name", "No. of genes", 
        "No. of significant pos. correlated genes",
        "No. of significant neg. correlated genes",
        "Median correlation strength")
    A[, 1] <- geneSetNames
    A[, 2] <- sizeGeneSets
    A[, 3] <- noOfSignCorrGenes
    A[, 4] <- noOfSignNegCorrGenes
    A[, 5] <- medianCorrelationStrength
    return(A)
}
