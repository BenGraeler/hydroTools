\name{hydroYears}
\alias{hydroYears}

\title{
Returns a list of hydrological years 
}
\description{
Generates a vector of characters identifying the hydrological years between start and end year. These strings will be understood by any time series from the \code{\link[=xts-package]{xts}} package.
}
\usage{
hydroYears(start, end, breaks = c("11-01", "10-31"))
}

\arguments{
  \item{start}{A single numeric specifying the start year.}
  \item{end}{A single numeric specifying the end year.}
  \item{breaks}{A vector of characters identifying month and day (MM-DD) of the annual breaks between hydrological periods (e.g. \code{breaks = c("11-01", "10-31")}).}
}

\value{A vector of characters identifying the temporal ranges of a series of hydrological periods. The format of the vector can be understood by a time series of the \code{\link[=xts-package]{xts}} package.}

\author{Benedikt Graeler}

\examples{
hydroYears(2013, 2014)

hydroYears(2011, 2013, breaks=c("11-01", "04-30", "05-01", "10-31"))
}