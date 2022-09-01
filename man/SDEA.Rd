\name{SDEA}
\alias{SDEA}

\title{
Super-Efficiency DEA
}
\description{
SDEA model to calculate the efficieny for each DMUs.
}
\usage{
SDEA(x=x, y=y, orientation = "input", rts = "crs", Cook = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ Inputs or resources used by each decision making unit.}
  \item{y}{ Outputs or resources used by each decision making unit.}
  \item{orientation}{ Orientation of the DEA model - primary emphasis on input-reduction input or output-augmentation output. Note the default is input orientation.}
  \item{rts}{ Returns to scale for the application, or industry studied. Note the default rts is crs. vrs Variable returns to scale. crs Constant returns to scale.}
  \item{Cook}{Used with variable returns to scale to address infeasibility in efficiency. Note the default is FALSE. Use TRUE when using vrs rts.}

}
\value{
\item{Input}{Input Values (x) passed to the model.}
\item{Output}{Output Values (y) passed to the model.}
\item{Orientation}{Orientation of the model.}
\item{RTS}{Returns to scale of the model.}
\item{Efficiency}{Efficiency of each DMU in the model.}
\item{Theta}{Used to calculate efficiency if the model is infeasiable. Note: Available only when Cook is set to TRUE.}
\item{Beta}{Used to calculate efficiency if the model is infeasiable.Note: Available only when Cook is set to TRUE.}
\item{Lambda}{Lambdas per DMU in the model.}
\item{StatusData}{Returns the status of the LP model.}
}
\references{
W.D. Cook, L. Liang, Y. Zha and J.Zhu (2009) A Modified Super-Efficiency DEA Model for Infeasibility, The Journal of the Operational Research Society Vol. 60, No. 2 (Feb., 2009), pp. 276-281.
}

\examples{
x <-data.frame(matrix(c(12, 26, 16, 60 ),ncol=1))
rownames(x) <- c('a','b','c','d')
y <- data.frame(matrix(c(6, 8, 9, 15 ),ncol=1))
rownames(y) <- c('a','b','c','d')


result <- SDEA(x=x, y=y, orientation = "input", rts = "crs", Cook = FALSE)
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
\keyword{ Two-Phase Approach }
\keyword{ Super Efficiency }
\keyword{ SDEA }
\keyword{ Malmquist Productivity Index }
\keyword{ MPI }% __ONLY ONE__ keyword per line

