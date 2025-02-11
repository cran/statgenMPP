## ----setup, include = FALSE-------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.dim = c(7, 4)
)
library(statgenMPP)
op <- options(width = 90, 
              digits = 3)

## ----simMarkerDat, results='asis', echo=FALSE-------------------------------------------
simMrkDat <- read.delim(system.file("extdata/multipop", "AxB.txt", package = "statgenMPP"))
knitr::kable(simMrkDat[1:4, 1:5])

## ----simPhenoDat, results='asis', echo=FALSE--------------------------------------------
simPhenoDat <- read.delim(system.file("extdata/multipop", "AxBxCpheno.txt", package = "statgenMPP"))
knitr::kable(head(simPhenoDat))

## ----simIBD-----------------------------------------------------------------------------
## Specify files containing markers.
# One file for each of the two crosses.
markerFiles <- c(system.file("extdata/multipop", "AxB.txt", 
                             package = "statgenMPP"),
                 system.file("extdata/multipop", "AxC.txt", 
                             package = "statgenMPP"))

## Specify file containing map.
# Both crosses use the same map file.
mapFile <- system.file("extdata/multipop", "mapfile.txt",
                       package = "statgenMPP")

## Read phenotypic data
phenoDat <- read.delim(system.file("extdata/multipop", "AxBxCpheno.txt",
                                   package = "statgenMPP"))
# Check contents.
head(phenoDat)

## Perform IBD calculations. 
ABCMPP <- calcIBDMPP(crossNames = c("AxB", "AxC"),  
                     markerFiles = markerFiles,
                     pheno = phenoDat,
                     popType = "F4DH",
                     mapFile = mapFile,
                     evalDist = 5)

## ----sumABCMPP--------------------------------------------------------------------------
## Print summary
summary(ABCMPP)

## ----plotPABCMPP------------------------------------------------------------------------
## Plot structure of the pedigree.
plot(ABCMPP, plotType = "pedigree")

## ----plotGABCMPP------------------------------------------------------------------------
## Plot genetic map.
# Highlight marker on chromosome 3 at position 40.
plot(ABCMPP, plotType = "genMap", highlight = "EXT_3_40")

## ----plotSGABCMPP, fig.height=10--------------------------------------------------------
## Plot IBD probabilities for genotype AxB0001.
plot(ABCMPP, plotType = "singleGeno", genotype = "AxB0001")

## ----ABCSQM-----------------------------------------------------------------------------
## Perform Single-QTL Mapping.
ABCSQM <- selQTLMPP(MPPobj = ABCMPP,
                    trait = "yield",
                    maxCofactors = 0)

## ----QPABCSQM---------------------------------------------------------------------------
## Plot QTL Profile for ABC SQM.
plot(ABCSQM, plotType = "QTLProfile")

## ----ABCMQM-----------------------------------------------------------------------------
## Perform Multi-QTL Mapping.
ABCMQM <- selQTLMPP(MPPobj = ABCMPP,
                    trait = "yield",
                    threshold = 3)

## ----ABCMQM_kin, eval=FALSE-------------------------------------------------------------
# ## Perform Multi-QTL Mapping.
# # Compute kinship matrices.
# ABCMQM_kin <- selQTLMPP(MPPobj = ABCMPP,
#                         trait = "yield",
#                         threshold = 3,
#                         computeKin = TRUE)

## ----plotQRABCMQM-----------------------------------------------------------------------
## Plot QTL Profile for ABC MQM.
plot(ABCMQM, plotType = "QTLRegion")

## ----plotQPABCMQM-----------------------------------------------------------------------
## Plot QTL Profile for ABC MQM.
plot(ABCMQM, plotType = "QTLProfile")

## ----plotPEABCMQM-----------------------------------------------------------------------
## Plot QTL Profile for ABC MQM.
plot(ABCMQM, plotType = "parEffs")

## ----plotQPEABCMQM----------------------------------------------------------------------
## Plot QTL Profile for ABC MQM.
plot(ABCMQM, plotType = "QTLProfileExt")

## ----plotCIABCMQM-----------------------------------------------------------------------
## Plot confidence intervals for parental effects for ABC MQM.
plot(ABCMQM, plotType = "parCIs")

## ----sumABCMQM--------------------------------------------------------------------------
## Print summary
summary(ABCMQM)

## ----extractABCRes----------------------------------------------------------------------
## Extract results of QTL mapping.
ABCMQMres <- ABCMQM$GWAResult$yield
head(ABCMQMres[, 1:8])

## ----extractABCQTL----------------------------------------------------------------------
## Extract QTLs and markers within QTL windows.
ABCMQMQTL <- ABCMQM$signSnp$yield
head(ABCMQMQTL[, c(2:8, 10)])

## ----maizeIBD---------------------------------------------------------------------------
## Define names of crosses.
crosses <- paste0("F353x", c("B73", "D06", "D09", "EC169", "F252", "F618", 
                             "Mo17", "UH250", "UH304", "W117"))
head(crosses)

## Specify files containing crosses.
## Extract them in a temporary directory.
tempDir <- tempdir()
crossFiles <- unzip(system.file("extdata/maize/maize.zip", package = "statgenMPP"), 
                    files = paste0(crosses, ".txt"), exdir = tempDir)

## Specify file containing map.
mapFile <- unzip(system.file("extdata/maize/maize.zip", package = "statgenMPP"), 
                 files = "map.txt", exdir = tempDir)

## Read phenotypic data.
phenoFile <- unzip(system.file("extdata/maize/maize.zip", package = "statgenMPP"),
                   files = "EUmaizePheno.txt", exdir = tempDir)
phenoDat <- read.delim(phenoFile)
head(phenoDat[, 1:5])

## Perform IBD calculations. 
maizeMPP <- calcIBDMPP(crossNames = crosses, 
                       markerFiles = crossFiles,
                       pheno = phenoDat,
                       popType = "DH",
                       mapFile = mapFile,
                       evalDist = 5)

## ----sumMaizeIBD------------------------------------------------------------------------
## Print summary
summary(maizeMPP)

## ----plotPmaizeIBD----------------------------------------------------------------------
## Plot structure of the pedigree.
plot(maizeMPP, plotType = "pedigree")

## ----maizeSQM, eval=FALSE---------------------------------------------------------------
# ## Perform Single-QTL Mapping.
# maizeSQM <- selQTLMPP(MPPobj = maizeMPP,
#                       trait = "mean_DtSILK",
#                       maxCofactors = 0)

## ----QPmaizeSQM-------------------------------------------------------------------------
## Plot QTL Profile for maize SQM.
plot(maizeSQM, plotType = "QTLProfile")

## ----maizeMQM, eval=FALSE---------------------------------------------------------------
# ## Perform Multi-QTL Mapping.
# maizeMQM <- selQTLMPP(MPPobj = maizeMPP,
#                       trait = "mean_DtSILK",
#                       threshold = 5)

## ----plotQPEmaizeMQM--------------------------------------------------------------------
## Plot QTL Profile for maize MQM.
plot(maizeMQM, plotType = "QTLProfileExt")

## ----plotCIamizeMQM---------------------------------------------------------------------
## Plot confidence intervals for parental effects for maize MQM.
plot(maizeMQM, plotType = "parCIs")

## ----barleyIBD--------------------------------------------------------------------------
## Specify files containing RABBIT output.
## Extract in a temporary directory.
tempDir <- tempdir()
inFile <- unzip(system.file("extdata/barley/barley_magicReconstruct.zip", 
                            package = "statgenMPP"), exdir = tempDir)

## Specify pedigree file.
pedFile <- system.file("extdata/barley/barley_pedInfo.csv",
                       package = "statgenMPP")

## Read phenotypic data.
data("barleyPheno")

## read RABBIT output. 
barleyMPP <- readRABBITMPP(infile = inFile,
                           pedFile = pedFile,
                           pheno = barleyPheno)

## ----sumPbarleyIBD----------------------------------------------------------------------
## Summary.
summary(barleyMPP)

## ----barleyMQM, eval=FALSE--------------------------------------------------------------
# ## Perform Multi-QTL Mapping with threshold 4.
# barleyMQM <- selQTLMPP(MPPobj = barleyMPP,
#                        trait = "Awn_length",
#                        threshold = 4)

## ----QPbarleyMQM16----------------------------------------------------------------------
## Plot QTL Profile for barley MQM - chromosome 1-6.
plot(barleyMQM, plotType = "QTLProfileExt", chr = 1:6)

## ----QPbarleyMQM7-----------------------------------------------------------------------
## Plot QTL Profile for barley MQM - chromosome 7.
plot(barleyMQM, plotType = "QTLProfileExt", chr = 7)

## ----plotCIbarleyMQM--------------------------------------------------------------------
## Plot confidence intervals for parental effects for maize MQM.
plot(barleyMQM, plotType = "parCIs")

## ----sumBarleyMQM-----------------------------------------------------------------------
## Summary.
summary(barleyMQM)

## ----ABCMQMPar, eval = FALSE------------------------------------------------------------
# ## Register parallel back-end with 2 cores.
# doParallel::registerDoParallel(cores = 2)
# 
# ## Perform Multi-QTL Mapping.
# ABCMQM_Par <- selQTLMPP(MPPobj = ABCMPP,
#                         trait = "yield",
#                         threshold = 3,
#                         parallel = TRUE)

## ----winddown, include = FALSE------------------------------------------------
options(op)

