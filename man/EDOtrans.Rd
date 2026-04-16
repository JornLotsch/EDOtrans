\name{EDOtrans}
\alias{EDOtrans}
\title{Euclidean distance-optimized data transformation}
\description{
  The package provides the necessary functions for performing the EDO data transformation.
}
\usage{
EDOtrans(Data, Cls, PlotIt = FALSE, FitAlg = "normalmixEM", Criterion = "LR",
                     MaxModes = 8, MaxCores = getOption("mc.cores", 2L), Seed = "simple")
}
\arguments{
 \item{Data}{the data as a vector.}
 \item{Cls}{the class information, if any, as a vector of similar length as instances in the data.}
 \item{PlotIt}{whether to plot the fit directly.}
 \item{FitAlg}{which fit algorithm to use: "ClusterRGMM" = GMM from ClusterR,
  "densityMclust" from mclust,
  "DO" from DistributionOptimization (slow),
  "MCMC" = NMixMCMC from mixAK, or
  "normalmixEM" from mixtools.}
 \item{Criterion}{which criterion should be used to establish the number of modes from the best GMM fit:
 "AIC", "BIC", "FM", "GAP", "LR" (likelihood ratio test), "NbClust" (from NbClust), "SI" (Silverman).}
 \item{MaxModes}{for automated GMM assessment: the maximum number of modes to be tried.}
 \item{MaxCores}{for automated GMM assessment: the maximum number of processor cores used under Unix.}
 \item{Seed}{Seed value for reproducibility. Options: \code{"auto"} (complex seed recovery, slow), \code{"simple"}
      (fast reproducible seed, default), or integer (exact control, recommended). Use integers for systematic testing with different seeds,
      \code{"simple"} for general use, or \code{"auto"} when you need to maintain exact RNG state continuity.}
}
\value{
  Returns a list of transformed data and class assignments.
 \item{DataEDO}{the EDO transformed data.}
 \item{EDOfactor}{the factor by which each data value has been divided.}
 \item{Cls}{the class information for each data instance.}
 }
\references{
Ultsch A, Lotsch J. Euclidean distance-optimized data transformation for cluster analysis in biomedical data (EDOtrans).
BMC Bioinformatics. 2022 Jun 16;23(1):233, https://doi.org/10.1186/s12859-022-04769-w
}
\author{
  Jorn Lotsch and Alfred Ultsch
}
\examples{
## example 1
data(iris)
IrisEDOdata <- EDOtrans(Data = as.vector(iris[,1]), Cls = iris$Species, Seed = 42)

## example 2
data(iris)
IrisEDOdata <-
  EDOtrans(Data = as.vector(iris[,1]), Cls = as.integer(iris$Species), MaxModes = 1, Seed = 42)

## example 3
data(iris)
set.seed(42)
IrisEDOdata <-
  EDOtrans(Data = as.vector(iris[,1]), Cls = iris$Species)
}

\keyword{EDOtrans}
\keyword{data transformation}
