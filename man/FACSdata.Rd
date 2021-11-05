\name{FACSdata}
\alias{FACSdata}

\docType{data}
\title{Example data of hematologic marker expression.}

\description{Data set of 4 flow cytometry-based lymphoma makers from 1559 cells from healthy subjects (class 1) and 1441 cells from lymphoma patients (class 2).}

\usage{data("FACSdata")}

\details{
Size 3000 x 4 , stored in \code{FACSdata$[FS,CDa,CDb,CDd]}
Original classes 2, stored in \code{FACSdata$Cls}
}

\examples{
data(FACSdata)
str(FACSdata)
}
