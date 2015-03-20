####################################################################
## Author: Bjarne Johannessen, Anita Sveen and Rolf I. Skotheim
## Maintainer: Bjarne Johannessen <bjarnej@ifi.uio.no>
## License: Artistic 2.0
## Part of the TIN package
####################################################################

## Function for creating a scatterPlot that compares the amount of splicing
## factor genes for which expression levels are significant positively and
## negatively correlated with the total relative amounts of aberrant exon 
## usage per sample.



posNegCorrPlot <- function(fileName, tra, geneSummaries, 
splicingFactors, randomGeneSets, traPermutations)
{

    fmat2 <- str_sub(fileName, -2, -1)
    fmat3 <- str_sub(fileName, -3, -1)
    fmat4 <- str_sub(fileName, -4, -1)

    spliceGeneSum <- 
        geneSummaries[rownames(geneSummaries) %in% splicingFactors[, 1], ]

    x <- vector('numeric', traPermutations)
    y <- vector('numeric', traPermutations)
    for (i in seq_len(traPermutations)) {
        perm <- sample(tra)
        corPerm <- corAndPvalue(t(spliceGeneSum[, seq_along(perm)]), 
            perm, alternative="two.sided")
        m <- matrix(c(corPerm$cor, corPerm$p), ncol = 2)
        colnames(m) <- c("cor", "p")
        m.df <- data.frame(m)
        s <- m.df[m.df$p < 0.05, ]
        x[i] <- length(s$cor[s$cor < 0])
        y[i] <- length(s$cor[s$cor > 0])
    }

    if (fmat3 == 'png') {
        png(fileName)
    } else if (fmat3 == 'jpg' || fmat4 == 'jpeg') {
        jpeg(fileName)
    } else if (fmat3 == 'eps' || fmat2 == 'ps') {
        postscript(fileName)
    } else if (fmat3 == 'pdf') {
        pdf(fileName)
    } else if (fmat3 == 'bmp') {
        bmp(fileName)
    }

    xx <- vector('numeric', randomGeneSets)
    yy <- vector('numeric', randomGeneSets)

    for (i in 1:randomGeneSets) {
        randGeneSummaries <- geneSummaries[sample(nrow(geneSummaries), 
            length(splicingFactors[, 1])), ]
        corPerm <- corAndPvalue(t(randGeneSummaries[, seq_along(tra)]), 
            tra, alternative="two.sided")
        g <- matrix(c(corPerm$cor, corPerm$p), ncol = 2)
        colnames(g) <- c("cor", "p")
        g.df <- data.frame(g)
        gs <- g.df[g.df$p < 0.05, ]
        xx[i] <- length(gs$cor[gs$cor < 0])
        yy[i] <- length(gs$cor[gs$cor > 0])
    }

    plot(x, y, pch = 20, col = "#2B3175", 
        xlim = c(0, length(splicingFactors[, 1])), 
        ylim = c(0, round(max(y, yy) + 25, digits = -1)), frame.plot = FALSE, 
        ylab = "No. of significant positively correlated genes", 
        xlab = "No. of significant negatively correlated genes", 
        axes = FALSE, cex = 2.0)
    axis(1, pos = 0, at = seq(0, length(splicingFactors[, 1]), by = 50), 
        labels = seq(0, length(splicingFactors[, 1]), by = 50))

    if (max(y, yy) > 50) {
        axis(2, pos = 0, at = seq(0, round(max(y, yy) + 50, digits = -1), 
            by = 50), labels = seq(0, round(max(y, yy) + 50, digits = -1), 
            by = 50))
    } else {
        axis(2, pos = 0, at = seq(0, 50, by = 10), 
            labels = seq(0, 50, by = 10))
    }
    points(xx, yy, pch = 20, col = "#7292CB", 
        xlim = c(0, length(splicingFactors[, 1])), cex = 2.0)

    correlationOrig <- corAndPvalue(t(spliceGeneSum[, 1:length(tra)]), 
        tra, alternative = "two.sided")
    m <- matrix(c(correlationOrig$cor, correlationOrig$p), ncol=2)
    colnames(m) <- c("cor", "p")
    m.df <- data.frame(m)
    gs <- m.df[m.df$p < 0.05, ]

    points(length(gs$cor[gs$cor<0]), length(gs$cor[gs$cor>0]), 
        pch = 20, col = "#9D2B49", cex = 3)

    path<-getwd()
    cat("Plot was saved in ",paste(path,"/",fileName,sep=""),"\n")
    dev.off()

}
