\name{Departments_Of_Accounting}
\alias{Departments_Of_Accounting}
\docType{data}
\title{
Data: UK University Departments Of Accounting Efficiency data.
}
\description{
Evaluation the Efficiency of UK University Departments Of Accounting Efficiency.
}
\usage{Departments_Of_Accounting}
\format{
  A data frame containing data for 20 UK University Departments Of Accounting.
  \describe{
    \item{\code{Departments}}{a numeric vector}
    \item{\code{Undergraduates}}{a numeric vector}
    \item{\code{Research}}{a numeric vector}
    \item{\code{Taught}}{a numeric vector}
    \item{\code{Res.Co}}{a numeric vector}
    \item{\code{OtherRes}}{a numeric vector}
    \item{\code{OtherIncome}}{a numeric vector}
    \item{\code{Publications}}{a numeric vector}
    \item{\code{AcademicStaff}}{a numeric vector}
    \item{\code{Salaries}}{a numeric vector}
    \item{\code{OtherExp}}{a numeric vector}
  }
}
\source{Tomkins C and Green RH (1988) An experiment in the use of data envelopment analysis for evaluating the efficiency of UK university departments of accounting. Financial Accounting and Management, 4, 147-164.}
\references{
Tomkins C and Green RH (1988) An experiment in the use of data envelopment analysis for evaluating the efficiency of UK university departments of accounting. Financial Accounting and Management, 4, 147-164.
}
\examples{
data(Departments_Of_Accounting)
attach(Departments_Of_Accounting)

x <- data.frame(AcademicStaff)
rownames(x) <- Departments
colnames(x) <- colnames(Departments_Of_Accounting)[9]

y <- data.frame(Undergraduates, Research, Taught,(Res.Co + OtherRes + OtherIncome))
rownames(y) <- Departments
colnames(y)[1] <- colnames(Departments_Of_Accounting)[2]
colnames(y)[2] <- colnames(Departments_Of_Accounting)[3]
colnames(y)[3] <- colnames(Departments_Of_Accounting)[4]
colnames(y)[4] <- c("Total_Income")

detach(Departments_Of_Accounting)

result <- DeaMultiplierModel(x,y,"crs", "input")

}
\keyword{datasets}
