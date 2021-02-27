\name{GMMartificialData}
\alias{GMMartificialData}

\docType{data}
\title{Example data an artificial Gaussioan mixture.}

\description{Dataset of 3000 instances with 3 variables that are Gaussian mixtures and belong to classes Cls = 1, 2, or 3, with different means and standard deviations and equal weights of 0.7, 0.3, and 0.1, respectively.}

\usage{data("GMMartificialData")}

\details{
Size 3000 x 3, stored in \code{GMMartificialData$[Var1,Var2,Var3]}

Classes 3, stored in \code{GMMartificialData$Cls}
}

\examples{
data(GMMartificialData)
str(GMMartificialData)
}
