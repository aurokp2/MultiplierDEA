\name{Japanese_Companies}
\alias{Japanese_Companies}
\docType{data}
\title{
Data: Japanese Companies data.
}
\description{
Japanese companies data for Operating Efficiency.
}
\usage{data("Japanese_Companies")}
\format{
  A data frame with 0 observations on the following 2 variables.
  \describe{
    \item{\code{DMU}}{a numeric vector}
    \item{\code{Company}}{a character vector}
    \item{\code{Asset}}{a numeric vector}
    \item{\code{Equity}}{a numeric vector}
    \item{\code{Employee}}{a numeric vector}
    \item{\code{Revenue}}{a numeric vector}
  }
}
\source{
W.D. Cook, L. Liang, Y. Zha and J.Zhu (2009) A Modified Super-Efficiency DEA Model for Infeasibility, The Journal of the Operational Research Society Vol. 60, No. 2 (Feb., 2009), pp. 276-281.
}
\references{
W.D. Cook, L. Liang, Y. Zha and J.Zhu (2009) A Modified Super-Efficiency DEA Model for Infeasibility, The Journal of the Operational Research Society Vol. 60, No. 2 (Feb., 2009), pp. 276-281.
}
\examples{
data(Japanese_Companies)
attach(Japanese_Companies)
}
\keyword{datasets}
