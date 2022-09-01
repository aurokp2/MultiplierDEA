\name{MPI}
\alias{MPI}

\title{
Malmquist Productivity Index.
}
\description{
MPI model to calculate MPI, Technical change, Efficiency change and Scale efficiency change.
}
\usage{
MPI(Dataset = Dataset, DMU_colName = DMU_colName, IP_colNames = IP_colNames, OP_ColNames = OP_ColNames, Period_ColName = Period_ColName, Periods = Periods, rts = "crs", orientation = "input", scale = FALSE)
}

\arguments{
  \item{Dataset}{The data required for the model.}
  \item{DMU_colName}{Column name for the DMUs in the dataset.}
  \item{IP_colNames}{Column name(s) for all input data in the dataset.}
  \item{OP_ColNames}{Column name(s) for all output data in the dataset.}
  \item{Period_ColName}{Column name for the period number in the dataset.}
  \item{Periods}{Unique periods numbers in the dataset in ascending order.}
  \item{rts}{Returns to scale for the application, or industry studied. Note the default rts is crs. vrs Variable returns to scale. crs Constant returns to scale.}
  \item{orientation}{Orientation of the DEA model - primary emphasis on input-reduction input or output-augmentation output. Note the default is input orientation.}
  \item{scale}{Note default value is FALSE.}
}

\value{
\item{DMU}{DMUs}
\item{et1t1.crs}{The efficiencies for period 1 with reference technology from period 1 for crs returns to scale. Note: available if returns to scale is crs or scale is TRUE.}
\item{et2t2.crs}{The efficiencies for period 2 with reference technology from period 2 for crs returns to scale. Note: available if returns to scale is crs or scale is TRUE.}
\item{et1t2.crs}{The efficiencies for period 2 with reference technology from period 1 for crs returns to scale. Note: available if returns to scale is crs or scale is TRUE.}
\item{et2t1.crs}{The efficiencies for period 1 with reference technology from period 2 for crs returns to scale. Note: available if returns to scale is crs or scale is TRUE.}
\item{et1t1.vrs}{The efficiencies for period 1 with reference technology from period 1 for vrs returns to scale. Note: available if returns to scale is vrs.}
\item{et2t2.vrs}{The efficiencies for period 2 with reference technology from period 2 for vrs returns to scale. Note: available if returns to scale is vrs.}
\item{et1t2.vrs}{The efficiencies for period 2 with reference technology from period 1 for vrs returns to scale. Note: available if returns to scale is vrs.}
\item{et2t1.vrs}{The efficiencies for period 1 with reference technology from period 2 for vrs returns to scale. Note: available if returns to scale is vrs}
\item{sec1}{First componenet of the scale efficiency change. (et1t2.crs/et1t2.vrs)/(et1t1.crs/et1t1.vrs)}
\item{sec2}{Second componenet of the scale efficiency change.(et2t2.crs/et2t2.vrs)/(et2t1.crs/et2t1.vrs)}
\item{sec}{Scale efficiency change. (sec1 * sec2) ^ 0.5}
\item{tc1}{First component of technical change. For crs, (et1t2.crs/et2t2.crs) and (et1t2.vrs/et2t2.vrs) for vrs.}
\item{tc2}{Second component of technical change. For crs, (et1t1.crs/et2t1.crs) and (et1t1.vrs/et2t1.vrs) for vrs.}
\item{tc}{Technical change. (tc1 * tc2) ^ 0.5}
\item{tec or ptec}{Efficiency change. Note: tec for crs and ptec for vrs returns to scale.}
\item{m.crs or m.vrs}{Malmquist Productivity index for the DMUs and periods.}
\item{Year}{Time period underconsideration for MPI.}
}

\references{
Rolf, Fare; Grosskopf, Shawna; Norris, Mary and Zhang, Zhongyang (1994) Productivity Growth, Technical Progress, and Efficiency Change in Industrialized Countries, The American Economic Review Vol. 84, No. 1, pp. 66-83.
}

\examples{

da_f <- data.frame(x= c(11, 29, 31, 61, 13, 27, 17, 61), y=c(6, 8, 11, 16, 7, 9, 10, 16), d= c(1,2,3,4, 1,2,3,4), p=c(1,1,1,1,2,2,2,2))

mpi_r <- MPI(Dataset = da_f, DMU_colName = "d", IP_colNames = "x", OP_ColNames = "y", Period_ColName = "p", Periods = c(1,2),rts = "vrs", orientation = "input", scale = TRUE)

# Examine the MPI for DMUs
mpi_r$m.vrs

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
