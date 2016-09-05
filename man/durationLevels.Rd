\name{durationLevels}
\alias{durationLevels}

\title{
Aggregation of a time series to set of duration levels.
}
\description{
A given time series is aggregated to a set of duration levels by means of a moving window.
}
\usage{
durationLevels(series, timeRes = 5/60, durLevels = c(1, 3, 6, 12, 24), 
               align="left", endrule = "NA", ..., NaN.rm = FALSE)
}

\arguments{
  \item{series}{A numeric vector holding the data of a times series.}
  \item{timeRes}{A single numeric value defining the temporal resolution of the time series with respect to the \code{durLevels} (e.g.: \code{timeRes = 5/60} if a 5-minute time series is given and \code{durLevels} are given in hours). }
  \item{durLevels}{A numeric vector of all the different duration levels that shall be used for the aggregation.}
  \item{align}{The alignment of the moving window, see \code{\link{runmean}}.}
  \item{endrule}{The endrule of the moving window, see \code{\link{runmean}}.}
  \item{...}{Arguments passed on to \code{\link{runmean}}.}
  \item{NaN.rm}{Shall \code{'NaN'} values be removed, i.e. set to \code{'NA'}?}
}
\details{
The moving window is realized using the function \code{\link{runmean}}, but the mean is multiplied by the temporal width to obtain the sum.  
}
\value{
A \code{\link{data.frame}} is returned where each column corresponds to a different aggregation level based on the argument \code{durLevels}. 
}
\author{
Benedikt Graeler}

\seealso{
\code{\link{runmean}}
}
\examples{
durLev <- durationLevels(runif(12*24*7))
dim(durLev)
colnames(durLev)
}
\keyword{aggregation}
\keyword{hydrology}