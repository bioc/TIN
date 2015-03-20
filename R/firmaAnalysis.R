####################################################################
## Author: Bjarne Johannessen, Anita Sveen and Rolf I. Skotheim
## Maintainer: Bjarne Johannessen <bjarnej@ifi.uio.no>
## License: Artistic 2.0
## Part of the TIN package
####################################################################

## Function for reading CEL files, and performing background correction, 
## normalization and alternative splicing analysis according to the FIRMA
## method.



firmaAnalysis <- function(useToyData=FALSE, aromaPath, dataSetName)
{

    if (useToyData) {
        sampleSetFirmaScores <- NULL
        data(sampleSetFirmaScores, envir=environment())
        fsScoresLog <- sampleSetFirmaScores
    } else {
        if (missing(aromaPath) && missing(dataSetName)) {
            stop("Please provide path to aroma.affymetrix root directory
                and name of data set as arguments",
                call.=FALSE)
        } else if (missing(aromaPath)) {
            stop("Please provide path to aroma.affymetrix root directory",
                call.=FALSE)
        } else if (missing(dataSetName)) {
            stop("Please provide name of data set",
                call.=FALSE)
        } else {
            oldPath<-getwd()
            setwd(aromaPath)
            chiptype <- "HuEx-1_0-st-v2"
            cdf <- AffymetrixCdfFile$byChipType(
                "HuEx-1_0-st-v2,coreR3,A20071112,EP")
            cs <- AffymetrixCelSet$byName(dataSetName, cdf=cdf)

            bc <- RmaBackgroundCorrection(cs)
            csBC <- process.RmaBackgroundCorrection(bc)

            qn <- QuantileNormalization(csBC, typesToUpdate="pm")
            csN <- process.QuantileNormalization(qn)

            plmTr <- ExonRmaPlm(csN, mergeGroups=TRUE)
            fit.ProbeLevelModel(plmTr)

            firma <- FirmaModel(plmTr)
            fit.FirmaModel(firma)
            fs <- getFirmaScores(firma)
            fsScores <- extractDataFrame(fs, units=NULL, addNames=TRUE)
            fsScoresLog <- log2(fsScores[, 6:length(fsScores)])
            setwd(oldPath)
        }
    }

    return(fsScoresLog)
}
