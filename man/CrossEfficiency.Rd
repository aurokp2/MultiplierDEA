\name{CrossEfficiency}
\alias{CrossEfficiency}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cross Efficiency Model
}
\description{
Cross Efficiency uses DEA to do peer evaluation of DMUs. Single-phase cross efficiency approach.
}
\usage{
CrossEfficiency(x = x, y = y, rts = "crs", orientation = "input", weightRestriction)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{x}{ Inputs or resources used by each decision making unit.}
  \item{y}{ Outputs or resources used by each decision making unit.}
  \item{rts}{ Returns to scale for the application, or industry studied. Note the default rts is crs. vrs Variable returns to scale. crs Constant returns to scale.}
  \item{orientation}{ Orientation of the DEA model - primary emphasis on input-reduction input or output-augmentation output. Note that unlike the DEA functions, the default is input orientation.}
  \item{weightRestriction}{ Weight restriction for the model. Optional parameter.}
}
\value{
The function returns a number of values per DMU.
\item{$ceva_matrix }{Returns the cross efficiency matrix. Row is the Rating DMU and Column is the Rated DMU.}
\item{$ce_ave }{Returns the cross efficiency score for the DMU.}
\item{$ceva_max }{Returns the maximum cross efficiency score for the DMU.}
\item{$ceva_min }{Returns the minimum cross efficiency score for the DMU.}
\item{$vx }{Input weights from the model.}
\item{$uy }{Output weights from the model.}
\item{$Model_Status }{Returns the status of the LP model.}
}
\note{
ceva_matrix - cross-evaluation matrix.
ceva_max - cross-evaluation maximum.
ceva_min - cross-evaluation minimum.
ce_ave - cross-efficiency scores.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
# Example from Kenneth R. Baker: Optimization Modeling with Spreadsheets, Third Edition, p. 176, John Wiley & Sons, Inc.

dmu <- c("A", "B", "C", "D", "E", "F")
x <- data.frame(c(150,400,320,520,350,320),c(0.2,0.7,1.2,2.0,1.2,0.7))
rownames(x) <- dmu
colnames(x)[1] <- c("StartHours")
colnames(x)[2] <- c("Supplies")

y <- data.frame(c(14,14,42,28,19,14),c(3.5,21,10.5,42,25,15))
rownames(y) <- dmu
colnames(y)[1] <- c("Reimbursed")
colnames(y)[2] <- c("Private")

# Calculate the efficiency score
result <- CrossEfficiency(x,y,"crs", "input")
# Examine the cross efficiency score for DMUs
print(result$ce_ave)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ DEA }
\keyword{ Multiplier DEA Model }
\keyword{ Cross Efficiency }
\keyword{ Single-Phase Approach }
\keyword{ Weight restrictions }
\keyword{ Assurance Region }
\keyword{ Malevolent }
\keyword{ Benevolent }
\keyword{ Two-Phase Approach }% __ONLY ONE__ keyword per line
