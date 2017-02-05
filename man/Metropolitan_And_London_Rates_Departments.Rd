\name{Metropolitan_And_London_Rates_Departments}
\alias{Metropolitan_And_London_Rates_Departments}
\docType{data}
\title{
Data: Metropolitan and London rates departments data
}
\description{
Relative Efficiency Metropolitan and London rates departments.
}
\usage{Metropolitan_And_London_Rates_Departments}
\format{
  A data frame containing data for 62 rates department authority.
  \describe{
    \item{\code{Authority}}{a character vector}
    \item{\code{TotalCost}}{a numeric vector}
    \item{\code{Non-cnl}}{a numeric vector}
    \item{\code{Rate}}{a numeric vector}
    \item{\code{Summons}}{a numeric vector}
    \item{\code{NPV}}{a numeric vector}
  }
}
\source{Dyson RG and Thanassoulis E (1988) Reducing weight flexibility in Data Envelopment Analysis, Journal of the Operational Research Society, 39 (6), 563-576.}
\references{
Dyson RG and Thanassoulis E (1988) Reducing weight flexibility in Data Envelopment Analysis, Journal of the Operational Research Society, 39 (6), 563-576.
}
\examples{
data(Metropolitan_And_London_Rates_Departments)
attach(Metropolitan_And_London_Rates_Departments)

x <- data.frame(TotalCost)
rownames(x) <- Authority
colnames(x) <- colnames(Metropolitan_And_London_Rates_Departments)[2]


y <- data.frame(`Non-cnl`, Rate, Summons, NPV)
rownames(y) <- Authority
colnames(y) <- colnames(Metropolitan_And_London_Rates_Departments)[3:6]

detach(Metropolitan_And_London_Rates_Departments)

result <- DeaMultiplierModel(x,y,"crs", "input")
}
\keyword{datasets}
