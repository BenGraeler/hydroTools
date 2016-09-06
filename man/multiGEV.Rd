\name{MultiGEV}
\alias{multiGEV}
\alias{dMultiGEV}
\alias{pMultiGEV}
\alias{qMultiGEV}
\alias{rMultiGEV}

\title{
The \code{multiGEV} Distribution
}

\description{
These functions provide information about a convex combination of GEV distributions. \code{dMultiGEV} gives the density, \code{pMultiGEV} gives the distribution function, \code{qMultiGEV} gives the quantile function and \code{rMultiGEV} generates random deviates.
}

\usage{
dMultiGEV(x, dur, par, log = F)
pMultiGEV(q, dur, par)
qMultiGEV(p, dur, par)
rMultiGEV(n, dur, par)

%rMultiGEV(n, dur, par, startId = NA)
}

\arguments{
  \item{x, q}{A numeric vector of quantiles.}
  \item{p}{A numeric vector of probabilities.}
  \item{n}{For simulation, the number of random numbers to draw.}
  \item{dur}{The duration levels to be evaluated.}
  \item{par}{The parameter list as returned by \code{\link{fitMultiGEV}}.}
  \item{log}{Single logical; shall the density be returned on log-scale?}
  % \item{startId}{Single numeric; starting index of the duration level for the simulation of random numbers. The default \code{startId = NA} causes a random starting index for each set of random numbers.}
}

\details{
The parameters of a set of GEV distributions are combined in convex manner for different duration levels. The resulting distribution family provides a smooth and non-decreasing distribution over distribution levels and increasing qunatiles. This ensures that events with larger return period and/or larger duration are always larger. 

The simulation is based on an fully independent sample. Hence, non-increasing sets of random numbers might occur.

% The simulation starts at a duration level according to \code{startId}. Moving away from the start index, random numbers are drawn from a restricted distribution such that the set of random numbers across the duration levels will always be increasing. 
}

\value{A vector or matrix corresponding densities, quantiles, fractions or set of random numbers for each duration.}

\author{
Benedikt Graeler}

\seealso{
\code{\link{dgev}}
}
\examples{
# set up initial (1) and increment (2) parameters of GEV for location sclae and shape
multiGEVpar <- c(50,5,5, # location
                 10,.5,.5, # scale
                 0,.5,.5) # shape

# set the base duration levels 
attr(multiGEVpar, "gevDl") <- c(1,6,24)

# plot qunatiles for different return periods and duration levels
contour(c(2:100), 1:24, qMultiGEV(1/2:100, dur = c(1:24), multiGEVpar),
        log="x", xlab="Time [a]", ylab="duration level [h]", ylim=c(0,24))
}

\keyword{ distribution }