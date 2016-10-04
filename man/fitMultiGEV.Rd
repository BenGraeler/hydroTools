\name{fitMultiGEV}
\alias{fitMultiGEV}

\title{A fitting routine for the multiGEV}

\description{
A fitting routine for the multiGEV
}

\usage{
fitMultiGEV(annMax, durLevels, start, parFun, ..., inf.res = 1e6)
}

\arguments{
  \item{annMax}{A numeric matrix; holding the annual maxima (rows) per duration level (columns).}
  \item{durLevels}{A nuemric vector; the duration levels used during the estimation.}
  \item{start}{A numeric vector; starting values for the parameter optimisation (\code{\link{optim}}).}
  \item{parFun}{A function that links the durations to the parameters of the GEV.}
  \item{\dots}{additional parmeters passed on to \code{\link{optim}}.}
  \item{inf.res}{A single numeric; the return value of theoptimized function if the log-likelihood happens to be infinite.}
}
\value{A numeric vector with parameters of the link function.}
\author{
Benedikt Graeler}

\seealso{\code{\link{multiGEV}}}

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
