\name{exQt}
\alias{exQt}

\title{
Counting exeedances of a given quantile in a moving window approach
}
\description{
Counting exeedances of a given quantile in a moving window approach
}
\usage{
exQt(ts, probs = 0.5, qt, w.width = 10, w.align = "right", exceed = TRUE, relative = TRUE)
}

\arguments{
  \item{ts}{A time series object.}
  \item{probs}{A numeric vector of quantile probabilities to consider.}
  \item{qt}{A numeric vector of quantiles; preferred over \code{probs}, but if missing empirically generated based on \code{probs}.}
  \item{w.width}{A single numeric defining the width of the window.}
  \item{w.align}{A single character specifying the alignment of the window (one of \code{"left"}, \code{"center"} or \code{"right"})}
  \item{exceed}{A single logical; whether the quantiles should be exceeded (the default) or undershot.}
  \item{relative}{A single logical; whether the count shall be weighted by the window width.}
}

\value{An \code{\link{xts}} object holding the counts. The columns either denote the set of window width, or the different quantiles.}

\author{Benedikt Graeler}

\examples{
library(xts)
set.seed(123)
dummyTs <- xts(rlnorm(12*24*7), as.POSIXct("2016-01-01 00:00:00") + 1:2016*5*60)
dummyTs["2016-01-02 00:00/2016-01-02 05:15",] <- 0
dummyTs["2016-01-02 12:25/2016-01-02 18:15",] <- 0
dummyTs["2016-01-04 16:45/2016-01-05 18:15",] <- 0

plot(dummyTs)

plot(exQt(dummyTs, 0.9, w.width = 20))
}
