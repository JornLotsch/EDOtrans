\name{CreateGMM}
\alias{CreateGMM}
\title{Creation of a Gaussian mixture data set}
\description{
  The function creates a data set according to a Gaussian mixture.
}
\usage{
CreateGMM(Means, SDs, Weights, n = 1000, Prob = FALSE)
}
\arguments{
 \item{Means}{a list of means of the mixture modes.}
 \item{SDs}{a list of standard deviations of the mixture modes.}
 \item{Weights}{a list of weights of the mixture modes, has to add up to 1.}
 \item{n}{number of data points.}
 \item{Prob}{whether a probability matrix should be calcuated.}
}
\value{
  Returns a data frame with the results, containing.
 \item{Data}{the generated data.}
 \item{Prob}{optional: the generated probabilities.}
 \item{Cls}{the classes of the instances.}
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
CreateGMM(Means = c(0,5,10), SDs = c(2,3,2), Weights = c(.4,.3,.3))
}
\keyword{EDOtrans}
\keyword{CreateGMM}
\keyword{data visualization}
\keyword{Gaussian mixture}
