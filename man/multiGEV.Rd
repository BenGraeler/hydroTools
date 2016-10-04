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
These functions provide a combination of GEV distributions based on a functional relationship between duration level and GEV parameters. \code{dMultiGEV} gives the density, \code{pMultiGEV} gives the distribution function, \code{qMultiGEV} gives the quantile function and \code{rMultiGEV} generates random deviates.
}

\usage{
dMultiGEV(x, dur, par = c(0, 1, 0), 
                  parFun = function(dur, par) par,
                  log = FALSE, oneColMat=FALSE)
pMultiGEV(q, dur, par = c(0, 1, 0), 
                  parFun = function(dur, par) par,
                  oneColMat=FALSE)
qMultiGEV(p, dur, par = c(0, 1, 0), 
                  parFun = function(dur, par) par,
                  oneColMat=FALSE)
rMultiGEV(n, dur, par = c(0, 1, 0),
                  parFun = function(dur, par) par,
                  oneColMat=FALSE)
}

\arguments{
  \item{x, q}{A numeric vector of quantiles.}
  \item{p}{A numeric vector of probabilities.}
  \item{n}{For simulation, the number of random numbers to draw per duration level \code{dur}.}
  \item{dur}{The duration levels to be evaluated.}
  \item{par}{The parameters of the link function \code{parFun}.}
  \item{parFun}{The link function providing the three GEV parameters for a given duration level.}
  \item{log}{Single logical; shall the density be returned on log-scale?}
  \item{oneColMat}{If \code{dur} has only length one, shall a one column matrix instead of a numeric vector be returned? (default = FALSE)}
}

\details{
The parameters of a set of GEV distributions are combined in the link function \code{parFun} for different duration levels. The resulting distribution family provides a smooth and non-decreasing distribution over distribution levels and increasing qunatiles. This ensures that events with larger return period and/or larger duration are always larger. 

The simulation is based on an fully independent sample. Hence, non-increasing sets of random numbers might occur.
}

\value{A vector or matrix of densities, quantiles, fractions or set of random numbers. Each column of the matrix corresponds to one duration level.}

\author{
Benedikt Graeler}

\seealso{
\code{\link{dgev}}
}
\examples{
library(evd)

# "manualy" sample a set of annual maxima
annMaxSim <- cbind(rgev(250, 10, 1, 0),
                   rgev(250, 15, 1.2, 0.1),
                   rgev(250, 20, 1.4, 0.2),
                   rgev(250, 25, 1.6, 0.3))

# define the parameter function returning the three GEV parameters for a given duration x
funPar <- function(x, par) 
  c(par[1] + par[2] * x,
    max(par[3] + par[4] * x, 0.01),
    par[5] + par[6] * x)

# fit a multiGEV to the simulated data assuming the relationship from "funPar"
fitPar <- fitMultiGEV(annMaxSim, 1:4*2, start =  c(5, 2.5, 0.8, 0.1, -0.1, 0.05), parFun = funPar)

# plot input parameters and fited parameters
par(mfrow=c(1,3))
plot(1:4*2, c(10, 15, 20, 25))
abline(fitPar[1], fitPar[2])

plot(1:4*2, c(1, 1.2, 1.4, 1.6))
abline(fitPar[3], fitPar[4])

plot(1:4*2, c(0, 0.1, 0.2, 0.3))
abline(fitPar[5], fitPar[6])
}

\keyword{ distribution }