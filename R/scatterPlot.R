####################################################################
## Author: Bjarne Johannessen, Anita Sveen and Rolf I. Skotheim
## Maintainer: Bjarne Johannessen <bjarnej@ifi.uio.no>
## License: Artistic 2.0
## Part of the TIN package
####################################################################

## Function for making a scatterplot showing relative amounts of aberrant
## exon usage per sample.



scatterPlot <-
    function(fileName, permutations=TRUE, percentileHits, permPercentileHits)
{
    fmat2 <- str_sub(fileName, -2, -1)
    fmat3 <- str_sub(fileName, -3, -1)
    fmat4 <- str_sub(fileName, -4, -1)

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
    } else {
        stop("Please use png, jpg, eps, pdf, or bmp as file format",
            call.=FALSE)
    }
    lows <- percentileHits[[1]]
    highs <- percentileHits[[2]]
    perm_lows <- permPercentileHits[[1]]
    perm_highs <- permPercentileHits[[2]]



    plot(log2(highs/ave(highs)), log2(lows/ave(lows)),
        xlim=c(-3, 3), ylim=c(-3, 3),
        type = "n", frame.plot=FALSE,
        ann=FALSE, axes=FALSE)
    axis(1, pos=0, at = -3:3, labels = c("-3", "-2", "-1", "", "1", "2", "3"))
    axis(2, pos=0, at = -3:3, labels = c("-3", "-2", "-1", "", "1", "2", "3"))
    text(2.0,-0.9, "Relative amounts of\naberrant inclusion (log2)")
    mtext("Relative amounts of\naberrant skipping (log2)", NORTH<-3, at=0)

    for(i in 1:length(highs)){
        x<-log2(highs[i]/ave(highs))
        y<-log2(lows[i]/ave(lows))
        points(x, y, pch = 20, col = "#5B83C2", cex = 1.5)
    }
    if (permutations) {
        points(log2(perm_highs/ave(perm_highs)),
            log2(perm_lows/ave(perm_lows)), pch=20, col='#E4A524')
    }
    path<-getwd()
    cat("Plot was saved in ",paste(path,"/",fileName,sep=""),"\n")
    dev.off()

}
