\name{Bank_Branch_Operating_Efficiency}
\alias{Bank_Branch_Operating_Efficiency}
\docType{data}
\title{
Data: Bank Branch Operating Efficiency data
}
\description{
Bank Branch data for Operating Efficiency.
}
\usage{Bank_Branch_Operating_Efficiency}
\format{
  A data frame containing data for 17 Bank Branches.
  \describe{
    \item{\code{Branch_Code}}{a character vector}
    \item{\code{PH}}{a numeric vector}
    \item{\code{OE}}{a numeric vector}
    \item{\code{SQM}}{a numeric vector}
    \item{\code{A}}{a numeric vector}
    \item{\code{B}}{a numeric vector}
    \item{\code{C}}{a numeric vector}
    \item{\code{D}}{a numeric vector}
  }
}
\source{Giokas DI (1991) Bank branck operating efficiency: A comparative application of DEA and the Loglinear model, OMEGA International Journal of Management Science, 19 (6) 549-557.
}

\references{
Giokas DI (1991) Bank branck operating efficiency: A comparative application of DEA and the Loglinear model, OMEGA International Journal of Management Science, 19 (6) 549-557.
}
\examples{
data(Bank_Branch_Operating_Efficiency)
attach(Bank_Branch_Operating_Efficiency)
x <- data.frame(PH, OE, SQM)
rownames(x) <- Branch_Code
colnames(x) <- colnames(Bank_Branch_Operating_Efficiency)[2:4]
y <- data.frame(A, B, C, D)
rownames(y) <- Branch_Code
colnames(y) <- colnames(Bank_Branch_Operating_Efficiency)[5:8]
detach(Bank_Branch_Operating_Efficiency)
# For CRS
result_CRS <- DeaMultiplierModel(x,y,"crs", "input")
# For VRS
result_VRS <- DeaMultiplierModel(x,y,"crs", "input")
}
\keyword{datasets}
