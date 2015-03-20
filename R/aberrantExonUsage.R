####################################################################
## Author: Bjarne Johannessen, Anita Sveen and Rolf I. Skotheim
## Maintainer: Bjarne Johannessen <bjarnej@ifi.uio.no>
## License: Artistic 2.0
## Part of the TIN package
####################################################################

## Function to calculate relative aberrant exon usage amounts per sample


aberrantExonUsage <- function(percentile, fs)
{

    quantile_low <- quantile(fs, probs=percentile/100, na.rm=TRUE)
    quantile_high <- quantile(fs, probs=1-(percentile/100), na.rm=TRUE)

    quantiles <- vector('list', 2)
    quantiles[[1]] <- quantile_low
    quantiles[[2]] <- quantile_high
    assign("quantiles", quantiles, envir = .GlobalEnv)

    llows <- vector('numeric', ncol(fs))
    lhighs <- vector('numeric', ncol(fs))

    for (i in seq_len(ncol(fs))) {
        llows[i] <- length(fs[fs[, i]<quantile_low, ][, i])
        lhighs[i] <- length(fs[fs[, i]>quantile_high, ][, i])
    }

    aberrantExons <- vector(mode='list')
    aberrantExons$skipping <- llows
    aberrantExons$inclusion <- lhighs
    assign("aberrantExons", aberrantExons, envir = .GlobalEnv)

    return(log2((lhighs + llows)/(ave(llows) + ave(lhighs))))

}
