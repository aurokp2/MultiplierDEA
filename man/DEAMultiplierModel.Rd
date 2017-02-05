\name{DeaMultiplierModel}
\alias{DeaMultiplierModel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
DEA Multiplier Model
}
\description{
DEA multiplier model calculates the efficieny and reference sets for each DMUs.
}
\usage{
DeaMultiplierModel(x = x, y = y, rts = "crs", orientation = "input", weightRestriction)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ Inputs or resources used by each decision making unit.}
  \item{y}{ Outputs or resources used by each decision making unit.}
  \item{rts}{ Returns to scale for the application, or industry studied. Note the default rts is crs. vrs Variable returns to scale. crs Constant returns to scale. Available option: crs, vrs}
  \item{orientation}{ Orientation of the DEA model - primary emphasis on input-reduction or output-augmentation output. Note that unlike the DEA functions, the default is input orientation. Available option: input, output.}
  \item{weightRestriction}{Weight restriction for the model. Optional parameter.}
}
\value{
The function returns a number of values per DMU. The standardized efficiency (all inefficiencies are between 0 and 1, for input and output orientation).
Efficiency, and lambda values are returned.
%%  If it is a LIST, use
\item{$rts }{Returns to scale of the model.}
\item{$Orientation }{Orientation of the model.}
\item{$InputValues }{Input Values (x) passed to the model.}
\item{$OutputValues }{Output Values (y) passed to the model.}
\item{$Efficiency }{Efficiency of each DMU in the model.}
\item{$Lambda }{Lambdas per DMU in the model.}
\item{$HCU_Input }{HCU data for inputs.}
\item{$HCU_Output }{HCU data for outputs.}
\item{$vx }{Input weights from the model.}
\item{$uy }{Output weights from the model.}
\item{$Free_Weights }{Free weights from the model. Applies only to vrs returns-to-scale.}
\item{$Model_Status }{Returns the status of the LP model.}

%% ...
}


%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
#Example from Kenneth R. Baker: Optimization Modeling with Spreadsheets, Third Edition,p. 176,
#John Wiley and Sons, Inc.

dmu <- c("A", "B", "C", "D", "E", "F")

x <- data.frame(c(150,400,320,520,350,320),c(0.2,0.7,1.2,2.0,1.2,0.7))
rownames(x) <- dmu
colnames(x)[1] <- c("StartHours")
colnames(x)[2] <- c("Supplies")

y <- data.frame(c(14,14,42,28,19,14),c(3.5,21,10.5,42,25,15))
rownames(y) <- dmu
colnames(y)[1] <- c("Reimbursed")
colnames(y)[2] <- c("Private")

#Creating the weight restriction data frame with Upper bound
weightRestriction<-data.frame(lower = c(1), numerator = c("StartHours"),
denominator = c("Supplies"), upper = c(2))

#Creating the weight restriction data frame without Upper bound
weightRestriction<-data.frame(lower = c(1), numerator = c("StartHours"),
denominator = c("Supplies"))

#Creating the weight restriction data frame with Upper bound and Na, Inf or NaN
weightRestriction<-data.frame(lower = c(1,2), numerator = c("StartHours","Reimbursed"),
denominator = c("Supplies","Private"), upper = c(2,Inf))


# Calculate the efficiency score without weight Restriction
result <- DeaMultiplierModel(x,y,"crs", "input")
# Examine the efficiency score for DMUs
print(result$Efficiency)


# Calculate the efficiency score with weight Restriction
result <- DeaMultiplierModel(x,y,"crs", "input", weightRestriction)
# Examine the efficiency score for DMUs
print(result$Efficiency)
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

