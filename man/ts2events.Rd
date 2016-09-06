\name{ts2events}
\alias{ts2events}

\title{
Event extraction from time series
}
\description{
A set of events is extracted from a time series based on a set of limiting conditions.
}
\usage{
ts2events(ts, var, AE = 0, tc = 0.5 * AE^0.6, 
          gap = max(4 * 60, 2 * tc * 60), 
          thd = 0.01, tunit = "mins")
}

\arguments{
  \item{ts}{A time series inheriting the class \code{\link{xts}}.}
  \item{var}{A character selecting the variable of interest, if missing, the first column is used.}
  \item{AE}{A single numeric defining the AE parameter (the area of the catchment).}
  \item{tc}{A single numeric defining the tc parameter (concentration time of the catchment).}
  \item{gap}{A single numeric defining the temporal gab between events (dry period).}
  \item{thd}{A single numeric setting the minimal threshold of the observed variable (precipitation amount).}
  \item{tunit}{A character defining the temporal resolution as used in \code{\link{difftime}} (e.g. \code{"mins"}).}
}

\value{A data.frame with columns \code{time}, \code{endTime}, \code{volume}, \code{peakTime}, \code{peak} and \code{duration} of the events.}

\references{
Maniak, U. (2005): Hydrologie und Wasserwirtschaft: Eine Einfuehrung fuer Ingenieure. Springer-Verlag.

Harms, R.W. (1984): Synthese massgebender Abflussdaten fuer teilbebaute kleine Einzugsgebiete. Mitteilungen des Instituts fuer Wasserwirtschaft der Universitaet Hannover, Heft 55, 1984.
}

\author{Benedikt Graeler}

\examples{
library(xts)
set.seed(123)
dummyTs <- xts(rlnorm(12*24*7), as.POSIXct("2016-01-01 00:00:00") + 1:2016*5*60)
dummyTs["2016-01-02 00:00/2016-01-02 05:15",] <- 0
dummyTs["2016-01-02 12:25/2016-01-02 18:15",] <- 0
dummyTs["2016-01-04 16:45/2016-01-05 18:15",] <- 0

plot(dummyTs)

ts2events(dummyTs, AE = 1)
}