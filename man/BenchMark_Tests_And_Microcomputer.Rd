\name{BenchMark_Tests_And_Microcomputer}
\alias{BenchMark_Tests_And_Microcomputer}
\docType{data}
\title{
Data: Relationship between benchmark tests and Microcomputer price data
}
\description{
The Relationship between benchmark tests and Microcomputer price data.
}
\usage{BenchMark_Tests_And_Microcomputer}
\format{
  A data frame containing data for 22 Microcomputers.
  \describe{
    \item{\code{System}}{a character vector}
    \item{\code{Price}}{a numeric vector}
    \item{\code{MemorySize}}{a numeric vector}
    \item{\code{DiskCapacity}}{a numeric vector}
    \item{\code{CPU}}{a numeric vector}
    \item{\code{IO}}{a numeric vector}
    \item{\code{RL1}}{a numeric vector}
    \item{\code{RL2}}{a numeric vector}
    \item{\code{RL3}}{a numeric vector}
  }
}
\source{Sircar S. and Dave D (1986) Tbe relationship between benchmark tests and microcomputer price. Communications of the ACM, 29, 212-217.
}
\references{
Sircar S. and Dave D (1986) Tbe relationship between benchmark tests and microcomputer price. Communications of the ACM, 29, 212-217.
}
\examples{
data(BenchMark_Tests_And_Microcomputer)
attach(BenchMark_Tests_And_Microcomputer)

x <- BenchMark_Tests_And_Microcomputer

detach(BenchMark_Tests_And_Microcomputer)

}
\keyword{datasets}
