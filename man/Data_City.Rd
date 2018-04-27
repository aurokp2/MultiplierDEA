\name{Data_City}
\alias{Data_City}
\docType{data}
\title{
Data: City data
}
\description{
City data for Operating Efficiency.
}
\usage{data("Data_City")}
\format{
  A data frame containing data for 15 city observations
  \describe{
    \item{\code{DMU}}{a numeric vector}
    \item{\code{City}}{a character vector}
    \item{\code{Houseprice}}{a numeric vector}
    \item{\code{Rental}}{a numeric vector}
    \item{\code{Violent}}{a numeric vector}
    \item{\code{Income}}{a numeric vector}
    \item{\code{B.Degree}}{a numeric vector}
    \item{\code{Doctor}}{a numeric vector}
  }
}
\source{
W.D. Cook, L. Liang, Y. Zha and J.Zhu (2009) A Modified Super-Efficiency DEA Model for Infeasibility, The Journal of the Operational Research Society Vol. 60, No. 2 (Feb., 2009), pp. 276-281.
}
\references{
W.D. Cook, L. Liang, Y. Zha and J.Zhu (2009) A Modified Super-Efficiency DEA Model for Infeasibility, The Journal of the Operational Research Society Vol. 60, No. 2 (Feb., 2009), pp. 276-281.
}
\examples{
data(Data_City)
attach(Data_City)
}
\keyword{datasets}
