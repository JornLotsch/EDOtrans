\name{GMMassessment}
\alias{GMMassessment}
\title{Gaussian mixture assessment}
\description{
  The function fits a Gaussian mixture and determines the number of statistically supported modes.
}
\usage{
GMMassessment(Data, DO = FALSE, PlotIt = FALSE, KS = FALSE, Criterion = "AIC", MaxModes = 10, MaxCores = 28, Seed)
}
\arguments{
 \item{Data}{the data as a vector.}
 \item{DO}{use the genetic algorithm distribution optimization or EM based GMM fitting (better but slow).}
 \item{PlotIt}{whether to plot the fit.}
 \item{KS}{perform a Kolmogorov-Smirnov test of the fit versus orginal distribution.}
 \item{Criterion}{criterion for mode number selection, likelihood ratio test LR, AIC or BIC.}
 \item{MaxModes}{the maximum number of modes to be tried.}
 \item{MaxCores}{the maximum number of processor cores used under Unix.}
 \item{Seed}{seed parameter set internally.}
}
\value{
  Returns a list of Gaussian modes.
 \item{Cls}{the classes to which the cases are assigned according to the Gaussian mode membership.}
 \item{Means}{means of the Gaussian modes.}
 \item{SDs}{standard deviations of the Gaussian modes.}
 \item{Weights}{weights of the Gaussian modes.}
 \item{Boundaries}{Bayesian boundaries between the Gaussian modes.}
 \item{plot}{Plot of the obtained mixture.}
 \item{KS}{Results of the Kolmogorov-Smirnov test.}
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
GMMassessment(Data = iris$Sepal.Length, DO = F, PlotIt = T, Criterion = "LR", MaxModes = 5, Seed = 42)
}
\keyword{EDOtrans}
\keyword{GMMassessment}
\keyword{data visualization}
