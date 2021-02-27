\name{EDOtrans}
\alias{EDOtrans}
\title{Euclidean distance-optimized data transformation}
\description{
  The package provides the necessary functions for performing the EDO data transformation.
}
\usage{
EDOtrans(Data, Cls, Means, SDs, Weights, DO = FALSE, PlotGMM = FALSE)
}
\arguments{
 \item{Data}{the data as a vector.}
 \item{Cls}{the class information, if any, as a vector of similar length as instances in the data.}
 \item{Means}{a list of mean values for a Gaussian mixure.}
 \item{SDs}{a list of standard deviations for a Gaussian mixure.}
 \item{Weights}{a list of weights for a Gaussian mixure.}
 \item{DO}{whether to use the DistributionOptimization method for Gaussian mixture modeling.}
 \item{PlotGMM}{whether to plot the obtained Gaussian mixture model.}
}
\value{
  Returns a list of transformed data and class assigments.
 \item{DataEDO}{the EDO transfomed data.}
 \item{EDOfactor}{the factor by which each data point has been deivded.}
 \item{Cls}{the class information for each data instance.}
 }
\references{
 Lotsch, J. (2021): 
 EDOtrans -- an R Package for Euclidean distance-optimized data transformation.
}
\author{
  Jorn Lotsch and Alfred Ultsch
}
\examples{
## example 1
data(iris)
IrisEDOdata <- EDOtrans(as.vector(iris[,1]), as.integer(iris$Species), PlotGMM = TRUE) 
}
\keyword{EDOtrans}
\keyword{data transformation}
