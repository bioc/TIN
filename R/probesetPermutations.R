####################################################################
## Author: Bjarne Johannessen, Anita Sveen and Rolf I. Skotheim
## Maintainer: Bjarne Johannessen <bjarnej@ifi.uio.no>
## License: Artistic 2.0
## Part of the TIN package
####################################################################

## function for creating permutations of the samples at each probeset.




probesetPermutations <- function(fs, percentiles)
{
    perm<-as.data.frame(t(apply(fs,1,sample)))

    perm_lows <- rep(NA, ncol(fs))
    perm_highs <- rep(NA, ncol(fs))

    for (i in 1:ncol(perm)) {
        perm_lows[i] <- length(perm[perm[i]<percentiles[1], ][, i])
        perm_highs[i] <- length(perm[perm[i]>percentiles[2], ][, i])
    }

    perms <- vector('list', 2)
    perms[[1]] <- perm_lows
    perms[[2]] <- perm_highs
    return(perms)

}
