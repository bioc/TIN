####################################################################
## Author: Bjarne Johannessen, Anita Sveen and Rolf I. Skotheim
## Maintainer: Bjarne Johannessen <bjarnej@ifi.uio.no>
## License: Artistic 2.0
## Part of the TIN package
####################################################################

## Function to create plot from hierarchical clustering analysis of the
## samples, based on splicing factor expression levels.



clusterPlot <- function(geneSummaries, tra, distmethod, clustermethod, 
    fileName)
{

    gs.t <- dist(t(geneSummaries), method=distmethod)

    h<-max(gs.t)
    l<-min(gs.t)
    gs.t <- 100*(gs.t-l+5)/(h-l)

    h.gs.t <- hclust(gs.t, method=clustermethod)

    ## Generate color map
    nyotin <- sapply(tra, function(x) max(x, -1.0))
    nyotin <- sapply(nyotin, function(x) min(x, 0.99))
    nyotin.cmap <- makecmap(nyotin, n = 30, 
        colFn = colorRampPalette(c('#6391C1', '#000000', '#9B3669')))

    nyotin.mat <- data.frame("oTIN" = cmap(nyotin, nyotin.cmap))
    rownames(nyotin.mat) <- colnames(geneSummaries)

    par(mar = c(5,4,4,3)+0.1)  # make space for color keys

    fmat2 <- str_sub(fileName, -2, -1)
    fmat3 <- str_sub(fileName, -3, -1)
    fmat4 <- str_sub(fileName, -4, -1)

    if (fmat3 == 'png') {
        png(fileName, width=480+max((length(tra)-20)*7, 0))
    } else if (fmat3 == 'jpg' || fmat4 == 'jpeg') {
        jpeg(fileName)
    } else if (fmat3 == 'eps' || fmat2 == 'ps') {
        postscript(fileName)
    } else if (fmat3 == 'pdf') {
        pdf(fileName)
    } else if (fmat3 == 'bmp') {
        bmp(fileName)
    }

    dendromat(h.gs.t, nyotin.mat, main = 'Hierarchical clustering', axes=FALSE)
    hkey(nyotin.cmap, side=3, 
        title = 'Total relative amounts of aberrant exon usage', stretch = 0.2)

    path<-getwd()
    cat("Plot was saved in ",paste(path,"/",fileName,sep=""),"\n")
    dev.off()

}
