\name{Evaluation_Educational_Program}
\alias{Evaluation_Educational_Program}
\docType{data}
\title{
Data: Educational program data
}
\description{
Evaluation of Educational program.
}
\usage{Evaluation_Educational_Program}
\format{
  A data frame containing data for 22 educational programs.
  \describe{
    \item{\code{Program}}{a numeric vector}
    \item{\code{CCR_EFF}}{a numeric vector}
    \item{\code{Revenue_Generated}}{a numeric vector}
    \item{\code{Student_Employed}}{a numeric vector}
    \item{\code{Employer_Satisfaction}}{a numeric vector}
    \item{\code{Contact_Hours}}{a numeric vector}
    \item{\code{Number_of_FTE_Staff}}{a numeric vector}
    \item{\code{Facility_Allocation}}{a numeric vector}
    \item{\code{Expenditures}}{a numeric vector}
  }
}
\source{Bessent A, Bessent W, Cbames A, Cooper WW and Thorgood N (1983) Evaluation of educational program proposals by means of data envelopment analysis. Educational Administrative Quarterly, 19, 82-107.}

\references{
Bessent A, Bessent W, Cbames A, Cooper WW and Thorgood N (1983) Evaluation of educational program proposals by means of data envelopment analysis. Educational Administrative Quarterly, 19, 82-107.
}
\examples{
data(Evaluation_Educational_Program)
attach(Evaluation_Educational_Program)


x <- data.frame(Contact_Hours, Number_of_FTE_Staff, Facility_Allocation, Expenditures)
rownames(x) <- Program
colnames(x) <- colnames(Evaluation_Educational_Program)[6:9]


y <- data.frame(Revenue_Generated, Student_Employed, Employer_Satisfaction)
rownames(y) <- Program
colnames(y) <- colnames(Evaluation_Educational_Program)[3:5]

detach(Evaluation_Educational_Program)

result <- DeaMultiplierModel(x,y,"crs", "input")
}
\keyword{datasets}
