\name{CrossEfficiency}
\alias{CrossEfficiency}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cross Efficiency Model
}
\description{
Cross Efficiency uses DEA to do peer evaluation of DMUs.
}
\usage{
CrossEfficiency<-function (x = x, y = y, rts = "crs", orientation = "input")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{x}{ Inputs or resources used by each decision making unit.}
  \item{y}{ Outputs or products of each decision making unit.}
  \item{rts}{ Returns to scale for the application, or industry studied. Note the default is crs rts. vrs Variable returns to scale. crs Constant returns to scale.}
  \item{orientation}{ Orientation of the DEA model - primary emphasis on input-reduction input or output-augmentation output. Note that unlike the DEA functions, the default is input orientation.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The function returns a number of values per DMU.
\item{$ce_matrix }{Returns the cross efficiency matrix. Row is the Rating DMU and Column is the Rated DMU.}
\item{$ce_ave }{Returns the cross efficiency score for the DMU.}
\item{$ce_max }{Returns the maximum cross efficiency score for the DMU.}
\item{$ce_min }{Returns the minimum cross efficiency score for the DMU.}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
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
print($ce_ave)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ DEA }
\keyword{ Multiplier DEA Model }
\keyword{ Cross Efficiency }% __ONLY ONE__ keyword per line
