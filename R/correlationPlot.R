####################################################################
## Author: Bjarne Johannessen, Anita Sveen and Rolf I. Skotheim
## Maintainer: Bjarne Johannessen <bjarnej@ifi.uio.no>
## License: Artistic 2.0
## Part of the TIN package
####################################################################

## Function for creating a plot that visualizes the number of splicing
## factor genes with expression levels significantly correlated with the
## sample-wise total relative amounts of aberrant exon usage.




correlationPlot <- function(fileName, tra, geneSummaries, 
    splicingFactors, randomGeneSets, traPermutations)
{

    fmat2 <- str_sub(fileName, -2, -1)
    fmat3 <- str_sub(fileName, -3, -1)
    fmat4 <- str_sub(fileName, -4, -1)

    spliceGeneSummaries <- 
        geneSummaries[rownames(geneSummaries) %in% splicingFactors[, 1], ]

    L <- vector('numeric', traPermutations)
    for (i in seq_len(traPermutations)) {
        perm <- sample(tra)
        cp <- corAndPvalue(t(spliceGeneSummaries[, seq_along(perm)]), 
            perm, alternative = "two.sided")
        m <- matrix(c(cp$cor, cp$p), ncol = 2)
        colnames(m) <- c("cor", "p")
        m.df <- data.frame(m)
        s <- m.df[m.df$p < 0.05, ]
        L[i] <- dim(s)[1]
    }

    bins <- seq(0, dim(splicingFactors)[1], by = 2)

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

    B <- vector('numeric', randomGeneSets)
    for (i in seq_len(randomGeneSets)) {
        randGeneSummaries <- geneSummaries[sample(nrow(geneSummaries),
            length(splicingFactors[, 1])), ]
        cp <- corAndPvalue(t(randGeneSummaries[, seq_along(tra)]), 
            tra, alternative = "two.sided")
        g <- matrix(c(cp$cor, cp$p), ncol = 2)
        colnames(g)<-c("cor", "p")
        g.df <- data.frame(g)
        gs <- g.df[g.df$p < 0.05, ]
        B[i] <- dim(gs)[1]
    }

    randomGeneSetsDist <- NULL
    traPermutationsDist <- NULL
    assign("randomGeneSetsDist", B, envir = .GlobalEnv)
    assign("traPermutationsDist", L, envir = .GlobalEnv)

    binsB <- seq(0, dim(splicingFactors)[1], by = 2)

    lim <- max(B, L)
    h <- hist(B, col = "#7292CB", border = "#7292CB", breaks = seq(min(B) - 2, 
        max(B) + 2, by = 2), axes = FALSE, main = "", 
        xlab = "No. of significantly correlated genes", 
        ylab = "No. of gene sets/permutations", 
        xlim = c(0, dim(splicingFactors)[1]))
    hh <- hist(L, breaks = seq(min(L) - 2, max(L) + 2, by = 2), 
        col = "#2B3175", border = "#2B3175", axes = FALSE)

    ylim = c(0,max(h$counts,hh$counts)) 
    nylim <- max(max(h$counts), max(hh$counts))

    plot(h, xlim = c(0, dim(splicingFactors)[1]), ylim = ylim, 
        col = "#7292CB", border = "#7292CB", axes = FALSE, 
        xlab = "No. of significantly correlated genes", 
        ylab = "No. of gene sets/permutations", main = "")
    plot(hh, xlim = c(0, dim(splicingFactors)[1]), ylim = ylim, add = TRUE, 
        col = "#2B3175", border = "#2B3175", axes=FALSE, main = "")

    axis(1, at = seq(0, dim(splicingFactors)[1], by = 70), 
        labels = seq(0, dim(splicingFactors)[1], by = 70))
    axis(2, at = c(0, round(max(hh$counts), digits = -1)), 
        labels = c(0, round(max(hh$counts), digits = -1)), las = 1)

    correlationOrig <- corAndPvalue(t(spliceGeneSummaries[, seq_along(tra)]), 
        tra, alternative = "two.sided")
    m <- matrix(c(correlationOrig$cor, correlationOrig$p), ncol = 2)
    colnames(m) <- c("cor", "p")
    m.df <- data.frame(m)
    s <- m.df[m.df$p < 0.05, ]

    points(dim(s)[1], 0, pch = 20, col = "#9D2B49", cex = 5)
    path<-getwd()
    cat("Plot was saved in ",paste(path,"/",fileName,sep=""),"\n")
    dev.off()
}
