\name{Evaluations_Of_NonProfitOrganizations}
\alias{Evaluations_Of_NonProfitOrganizations}
\docType{data}
\title{
Data: Evaluation of Non-Profit organizations data
}
\description{
Evaluation of Non-Profit organizations efficiency.
}
\usage{Evaluations_Of_NonProfitOrganizations}
\format{
  A data frame containing data for 16 Non-Profit organizations.
  \describe{
    \item{\code{Hospital}}{a numeric vector}
    \item{\code{H0}}{a numeric vector}
    \item{\code{PercentOccupancy}}{a numeric vector}
    \item{\code{RevenuePerDay}}{a numeric vector}
    \item{\code{A/RTurnover}}{a numeric vector}
    \item{\code{CostPerDay}}{a numeric vector}
    \item{\code{LengthOfStay}}{a numeric vector}
    }
}
\source{Greenberg R and Nunamaker T (1987) A generalized multiple criteria model for control and evaluation of nonprofit organizations. Financial Accountability and Management, 3 (4), 331-342.}
\references{
Greenberg R and Nunamaker T (1987) A generalized multiple criteria model for control and evaluation of nonprofit organizations. Financial Accountability and Management, 3 (4), 331-342.
}
\examples{
data(Evaluations_Of_NonProfitOrganizations)
attach(Evaluations_Of_NonProfitOrganizations)

x <- Evaluations_Of_NonProfitOrganizations

detach(Evaluations_Of_NonProfitOrganizations)

}
\keyword{datasets}
