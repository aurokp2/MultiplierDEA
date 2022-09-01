\name{Mal_Ben}
\alias{Mal_Ben}
\alias{BenAndMal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Benevolent and Malevolent Model
}
\description{
Two-Phase Cross efficiency approach.
}
\usage{
Mal_Ben(x = x, y = y, rts ="crs", orientation = "input", phase = "mal",
weightRestriction, include = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{x}{ Inputs or resources used by each decision making unit.}
  \item{y}{ Outputs or resources used by each decision making unit.}
  \item{rts}{ Returns to scale for the application, or industry studied. Note the default rts is crs. vrs Variable returns to scale. crs Constant returns to scale. Available option: crs, vrs.}
  \item{orientation}{ Orientation of the DEA model - primary emphasis on input-reduction input or output-augmentation output. Note that unlike the DEA functions, the default is input orientation. Available option: input, output.}
  \item{weightRestriction}{Weight restriction for the model. Optional parameter.}
  \item{phase}{Second phase of the model. Malevolent or Benevolent. Note the default is mal.Available option: mal, ben.}
  \item{include}{In the second phase include evaluating DMU in the calculation. Default is TRUE. Available option: TRUE, FALSE.}
}
\value{
The function returns a number of values per DMU. The standardized efficiency (all inefficiencies are between 0 and 1, for input and output orientation)
Efficiency, and the lambda values, lambda, are returned.
%%  If it is a LIST, use
\item{$rts }{Returns to scale of the model.}
\item{$Orientation }{Orientation of the model.}
\item{$InputValues }{Input Values (x) passed to the model.}
\item{$OutputValues }{Output Values (y) passed to the model.}
\item{$Phase1_Efficiency }{Efficiency of each DMU in the model from Phase 1.}
\item{$Phase1_Lambda }{Lambdas per DMU in the model from Phase 1.}
\item{$Phase1_vx }{Input weights from the model from Phase 1.}
\item{$Phase1_uy }{Output weights from the model from Phase 1.}
\item{$Phase1_Free_Weights }{Free weights from the model from Phase 1. Applies only to vrs returns-to-scale.}
\item{$Phase1_Model_Status }{Returns the status of the phase two LP model.}
\item{$Phase2_Efficiency }{Efficiency of each DMU in the model from Phase 2.}
\item{$Phase2_Lambda }{Lambdas per DMU in the model from Phase 2.}
\item{$Phase2_vx }{Input weights from the model from Phase 2.}
\item{$Phase2_uy }{Output weights from the model from Phase 2.}
\item{$Phase2_Free_weights }{Free weights from the model from Phase 2. Applies only to vrs returns-to-scale.}
\item{$Phase2_Model_Status }{Returns the status of the phase two LP model. }
\item{$ceva_matrix }{Returns the cross efficiency matrix. Row is the Rating DMU and Column is the Rated DMU.}
\item{$ce_ave }{Returns the cross efficiency score for the DMU.}
\item{$ceva_max }{Returns the maximum cross efficiency score for the DMU.}
\item{$ceva_min }{Returns the minimum cross efficiency score for the DMU.}
}
\note{
ceva_matrix - cross-evaluation matrix.
ceva_max - cross-evaluation maximum.
ceva_min - cross-evaluation minimum.
ce_ave - cross-efficiency scores.
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

# Calculate the efficiency score
result <- Mal_Ben(x, y, rts = "crs",  orientation = "input",phase = "mal", include = TRUE)
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
\keyword{ Two-Phase Approach }
\keyword{ Super Efficiency }
\keyword{ SDEA }
\keyword{ Malmquist Productivity Index }
\keyword{ MPI }% __ONLY ONE__ keyword per line
